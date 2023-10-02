#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");


my $usage = "usage: $0 fusion_result_file_listing.dat\n\n";

my $fusion_result_file_listing = $ARGV[0] or die $usage;

my $fusion_prog_parser_lib_dir = "$FindBin::Bin/../../benchmarking/FusionProgParsers";

my %prog_type_to_file_parser = ( 

    'ctat-LR-fusion.v0.10.0' => 'CTAT_LR_parser',
    'ctat-LR-fusion.v0.10.2' => 'CTAT_LR_parser',
    'flairfusion' => 'FlairFusion_parser',
    'fusionseeker' => 'FusionSeeker_parser',
    'JAFFAL' => 'JAFFAL_parser',
    'LongGF' => 'LongGF_parser',
    'pbfusion_v0.3.1' => 'PBfusion_parser',
    );


foreach my $module (values %prog_type_to_file_parser) {
    my $module_path = "$fusion_prog_parser_lib_dir/$module.pm";

    require($module_path);

}


main: {

    # print header
    print join("\t", "sample", "prog", "fusion", "breakpoint", "num_reads") . "\n";
    
    open(my $fh, $fusion_result_file_listing) or die "Error, cannot open file $fusion_result_file_listing";
    while (<$fh>) {
        my $line = $_;
        #print $line;
        chomp;
                
        my ($prog_name, $dataset, $result_file) = split(/\t/);
        
        my $parser_module;

        if (exists $prog_type_to_file_parser{$prog_name}) {
            $parser_module = $prog_type_to_file_parser{$prog_name};
        }

        unless (defined $parser_module) {
            die "Error, no parser for prog [$prog_name] ";
        }
        
        my $parser_function = $parser_module . "::" . "parse_fusion_result_file";
        
        no strict 'refs';
        my @fusions = &$parser_function($result_file);

        #use Data::Dumper;
        #print(Dumper(\@fusions));
        
        unless (@fusions) {
            die "Error, no fusions reported for $line";
        }
        
        @fusions = reverse sort { $a->{num_reads} <=> $b->{num_reads} } @fusions;
        
        foreach my $fusion (@fusions) {

            my $fusion_name = join("--", $fusion->{geneA}, $fusion->{geneB});
            
            my $num_reads = $fusion->{num_reads};
            
            my $breakpoint = $fusion->{breakpoint};
            
            print join("\t", $dataset, $prog_name, $fusion_name, $breakpoint, $num_reads) . "\n";
        }
        
    }
    close $fh;


    exit(0);
}

