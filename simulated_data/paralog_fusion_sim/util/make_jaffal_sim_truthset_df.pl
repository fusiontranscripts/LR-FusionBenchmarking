#!/usr/bin/env perl

use strict;
use warnings;


my %read_breakpoint_info;
# must run in the LongReadFusionSimulation/ dir.
open (my $fh, "../fusion_breakpoint_summary.tsv") or die $!;
my $header = <$fh>;
chomp $header;
while (<$fh>) {
    my $line = $_;
    chomp $line;
    my @vals = split(/\t/, $line);
    my $descr_line = $vals[8];
    my @pts = split(/\s+/, $descr_line);
    my $acc = $pts[0];
    $read_breakpoint_info{$acc} = $line;
}
close $fh;


print "$header\tsample\tfusion_id\tnum_reads\n";

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
        my $breakpoint_info = $read_breakpoint_info{$fusion} or die "Error, no breakpoint info for $fusion";
        print join("\t", $breakpoint_info, $dataset, $fusion, $num_reads) . "\n";
    }
}

exit(0);
