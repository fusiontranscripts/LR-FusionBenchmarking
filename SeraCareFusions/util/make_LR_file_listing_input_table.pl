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

    if ($filename =~ /bc08/i) {
        $dataset = "MAS-seq-L1";
    }
    elsif ($filename =~ /bc09/) {
        $dataset = "MAS-seq-L2";
    }
    elsif ($filename =~ /m64020e_230609_061545|isoseq/i) {
        $dataset = "ISO-seq";
    }
    else {
        die "Error, cannot ascertain data set type for $dataset of $filename";
    }
    
        
    if (%restrict_progs && ! exists $restrict_progs{$progname}) {
        print STDERR "make_file_listing_input_table::  - skipping $filename, not in restricted list.\n";
        next;
    }
    
        
    print join("\t", $progname, $dataset, $filename) . "\n";
    
    #my $proper_progname = $converter{$prog};

     
    # keep original name
    #print join("\t", $sample_name, $prog, $filename) . "\n";
    
}
    



exit(0);


    
