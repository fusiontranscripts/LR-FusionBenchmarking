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


#    if ($progname eq "JAFFAL") {
#        if ($filename =~ /jaffa_results.csv/) {
#            print join("\t", "JAFFAL", "all", $filename) . "\n";
#        }
#        next;
#    }
    
#    if ($progname =~ /longgf/i && $dataset !~ /fa_longgf.logs$/) {
#        next;
#    }
    
#    if ($filename =~ /fusionseeker/) {
#        $dataset = join("/", $pts[-2], $pts[-1]);
#        $progname = $pts[-3];
#        if ($filename !~ /confident_genefusion.txt/) {
#            next;
#        }
#    }
    
    
    
    if (%restrict_progs && ! exists $restrict_progs{$progname}) {
        print STDERR "make_file_listing_input_table::  - skipping $filename, not in restricted list.\n";
        next;
    }
    
    my $seqtype;
    if ($dataset =~ /Pac/) {
        $seqtype = "Pac";
    }
    elsif ($dataset =~ /ONT/) {
        $seqtype = "ONT";
    }
    else {
        die "Error, cannot ascertain Pac or ONT from $dataset";
    }

    my $divergence_level;
    if ($dataset =~ /(\d+)err/) {
        $divergence_level = $1;
    }
    else {
        die "Error, cannot ascertain divergence level from $dataset";
    }
    
    
    $dataset = "${seqtype}_${divergence_level}";
    
    print join("\t", $progname, $dataset, $filename) . "\n";
    
    #my $proper_progname = $converter{$prog};

     
    # keep original name
    #print join("\t", $sample_name, $prog, $filename) . "\n";
    
}
    



exit(0);


    
