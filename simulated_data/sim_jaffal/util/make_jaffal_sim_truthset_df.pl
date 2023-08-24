#!/usr/bin/env perl

use strict;
use warnings;


print "sample\tfusion\tnum_reads\n";

foreach my $file (<*.truthset>) {
    
    my $seqtype;
    if ($file =~ /Pac/) {
        $seqtype = "Pac";
    }
    elsif ($file =~ /ONT/) {
        $seqtype = "ONT";
    }
    else {
        die "Error, cannot ascertain Pac or ONT from $file";
    }

    my $divergence_level;
    if ($file =~ /(\d+)err/) {
        $divergence_level = $1;
    }
    else {
        die "Error, cannot ascertain divergence level from $file";
    }
    
    
    my $dataset = "${seqtype}_${divergence_level}";
    
    

    open(my $fh, $file) or die $!;
    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my ($fusion, $num_reads) = split(/\t/);
        print join("\t", $dataset, $fusion, $num_reads) . "\n";
    }
}

exit(0);
