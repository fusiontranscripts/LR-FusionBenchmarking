#!/usr/bin/env perl

use strict;
use warnings;
use Carp;


my %restrict_progs;
if (@ARGV) {
    my $restrict_progs_file = $ARGV[0];
    open(my $fh, $restrict_progs_file) or die "Error, cannot open file: $restrict_progs_file";
    while(<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my $progname = $_;
        $restrict_progs{$progname} = 1;
    }
    close $fh;
}
    

## convert prog name tokens to names used in the data table.
my %converter = ( );


while (<STDIN>) {
    chomp;
    my $filename = $_;
    chomp $filename;
    
    unless (-f $filename) {
        #print STDERR "warning, $filename is not a file. Skipping...\n";
        next;
    }
    
    my @pts = split(/\//, $filename);
    my $dataset = $pts[-1];
    my $progname = $pts[-2];


    if ($filename =~ /fusionseeker/) {
        $dataset = join("/", $pts[-2], $pts[-1]);
        $progname = $pts[-3];
        if ($filename !~ /confident_genefusion.txt/) {
            next;
        }
    }
    
    
    
    if (%restrict_progs && ! exists $restrict_progs{$progname}) {
        print STDERR "make_file_listing_input_table::  - skipping $filename, not in restricted list.\n";
        next;
    }


    my $cov;
    my $sample;
    
    if ($filename =~ /ont_(cov\d+)_(sample\d+)/) {
        $cov = $1;
        $sample = $2;
    }
    else {
        die "Error, cannot extract cov and sample info from filename: $filename";
    }

    $dataset = "${cov}_${sample}";
    
    print join("\t", $progname, $dataset, $filename) . "\n";
        
}

exit(0);


    
