#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 consolidated_edgren_predictions.dat min_agree\n\n";

my $input_file = $ARGV[0] or die $usage;
my $min_agree = $ARGV[1] or die $usage;

main: {
    


    my %fusion_to_prog;
    my %orig_fusion_call;
    my %prognames;

    open(my $fh, $input_file) or die $!;
    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $score_class = $x[0];
        my $sample_name = $x[1];
        my $progname = $x[2];
        my $fusion_name = $x[9];
        $fusion_name = "$sample_name|$fusion_name";
        
        $prognames{$progname}++;
        
        if ($score_class eq "TP") {
            $fusion_to_prog{$fusion_name}->{$progname}++;
        }            
    }
    
    ## capture those fusions that meet the min prog criteria

    my %fusions_meet_min_prog_count;

    foreach my $fusion_name (keys %fusion_to_prog) {
        
        my $prog_count = scalar(keys %{$fusion_to_prog{$fusion_name}});
        
        #print "$fusion_name\t$orig_fusion_name\t$prog_count\n";

        if ($prog_count >= $min_agree) {
            $fusions_meet_min_prog_count{$fusion_name} = 1;
        }
    }
        

    ## generate report
    my @prognames = sort keys %prognames;

    print "prog\t" . join("\t", @prognames) . "\n";

    my @final_fusions = sort keys %fusions_meet_min_prog_count; 
    
    foreach my $fusion (@final_fusions) {
        
        my @vals = ($fusion);
        foreach my $progname (@prognames) {
            my $found = (exists $fusion_to_prog{$fusion}->{$progname}) ? 1 : 0;
            push (@vals, $found);
        }

        print join("\t", @vals) . "\n";
    }

    exit(0);
    
}
