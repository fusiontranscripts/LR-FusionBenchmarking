#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


=descr

Removing messy fusions defined as:
- gene is reported as part of a fusion by multiple programs across multiple samples.

=cut



my $usage = "usage: $0 fusion_preds.tsv progs_select.txt\n\n";

my $fusion_preds_file = $ARGV[0] or die $usage;
my $progs_select_file = $ARGV[1] or die $usage;

my %progs_select;
open(my $fh, $progs_select_file) or die $!;
while(<$fh>) {
    s/\s+//;
    unless (/\w/) { next; }
    if (/^\#/) {
        next;
    }
    $progs_select{$_} = 1;
}
close $fh;

open($fh, $fusion_preds_file) or die $!;
my $delim_parser = new DelimParser::Reader($fh, "\t");

my @column_headers = $delim_parser->get_column_headers();

my %gene_to_prog;

my @rows;
while(my $row = $delim_parser->get_row()) {
    
    push (@rows, $row);
        
    my $prog = $row->{prog};
    my $sample = $row->{sample};
    
    if ($progs_select{$prog}) {
        # include in scoring
        foreach my $gene (split(/--/, $row->{fusion})) {
            $gene_to_prog{$gene}->{$prog}->{$sample} = 1;
        }
    }
}



open (my $ofh_pass, ">$fusion_preds_file.pass") or die $!;
open(my $ofh_messy, ">$fusion_preds_file.messy_fail") or die $!;

my $pass_writer = new DelimParser::Writer($ofh_pass, "\t", \@column_headers);
my $messy_writer = new DelimParser::Writer($ofh_messy, "\t", [@column_headers, "messy_reason"]);

foreach my $row (@rows) {
    
    my @genes = split(/--/, $row->{fusion});
        
    if (my $messy_gene_reasons = &examine_for_messy_gene(@genes) ) {
        
        $row->{messy_reason} = $messy_gene_reasons;
        $messy_writer->write_row($row);
    }
    else {
        $pass_writer->write_row($row);
    }
    
}


exit(0);


####
sub examine_for_messy_gene {
    my @genes = @_;

    my @messy;
    
    foreach my $gene (@genes) {
        
        my @progs = keys %{$gene_to_prog{$gene}};
        if (scalar @progs > 1) {
            # multiple progws
            foreach my $prog (@progs) {
                my @samples = keys %{$gene_to_prog{$gene}->{$prog}};
                if (scalar(@samples) >= 3) {
                    # and multiple samples
                    push (@messy, "$gene^$prog^" . join(",", @samples));
                }
            }
        }
    }
    
    if (scalar(@messy) > 1) {
        # multiple progs, each mult samples
        
        return(join(";", @messy));
    }
    
    else {
        return;
    }
}
