#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 ref_annot_gtf.gene_spans\n\n";

my $gene_spans_file = $ARGV[0] or die $usage;

=gene_spans_file_format
ENSG00000225538.1       chr11   56082801        56083739        +       OR5BE1P unprocessed_pseudogene
ENSG00000237851.1       chr6    142788123       142794086       +       RP1-67K17.4     lincRNA
ENSG00000212855.5       chrY    9740584 9758476 +       TTTY2   lincRNA
ENSG00000257527.1       chr16   18411309        18411851        -       MIR3179-3       lincRNA
ENSG00000278179.1       chr5    37888249        37888350        -       AC008869.2      miRNA
ENSG00000130529.14      chr19   49157741        49211836        +       TRPM4   protein_coding
ENSG00000225193.5       chr15   90206572        90206967        +       RPS12P26        transcribed_processed_pseudogene
=cut


my %gene_symbol_counter;
open(my $fh, $gene_spans_file) or die $!;

# first count gene symbols:

while(<$fh>) {
    chomp;
    my ($ensg_id, $chrom, $lend, $rend, $orient, $gene_symbol, $gene_type) = split(/\t/);
    $gene_symbol_counter{$gene_symbol}++;
}

open($fh, $gene_spans_file) or die $!;
print join("\t", "gene_id",	"chr",	"lend",	"rend",	"file",	"primary_target") . "\n";
while(<$fh>) {
    chomp;
    my ($ensg_id, $chrom, $lend, $rend, $orient, $gene_symbol, $gene_type) = split(/\t/);
    # format want: 
    # gene_id chr     lend    rend    file    primary_target

    if ($gene_symbol_counter{$gene_symbol} > 3) { next; } # avoid SNORNAs and similar
    
    print join("\t", $gene_symbol, $chrom, $lend, $rend, "ref_annot_gtf.gene_spans", 1) . "\n";
}

    
