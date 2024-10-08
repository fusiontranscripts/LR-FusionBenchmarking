#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;

my $usage = "usage: $0 preds.collected\n\n";

my $preds_file = $ARGV[0] or die $usage;


main: {

    my %prognames;
    my %fusion_preds;

    open(my $fh, $preds_file) or die $!;
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    
    while (my $row = $delim_parser->get_row()) {
        
        my $sample_name = $row->{sample};
        my $prog = $row->{prog};
        my $fusion_name = uc $row->{fusion};

        my $num_reads = $row->{num_reads};
        
        $fusion_name = "$sample_name|$fusion_name";

        $prognames{$prog} = 1;
        
        $fusion_preds{$fusion_name}->{$prog} = $num_reads;
        
    }
    close $fh;


    ## output matrix
    my @prognames = sort keys %prognames;
    my @fusions = sort keys %fusion_preds;

    print "\t" . join("\t", @prognames) . "\n";

    foreach my $fusion (@fusions) {
        my @vals = ($fusion);
        foreach my $progname (@prognames) {
            my $val = $fusion_preds{$fusion}->{$progname} || 0;
            push (@vals, $val);
        }
        
        print join("\t", @vals) . "\n";
    }

    exit(0);
}





