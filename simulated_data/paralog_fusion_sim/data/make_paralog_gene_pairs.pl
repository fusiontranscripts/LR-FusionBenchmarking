#!/usr/bin/env perl

use strict;
use warnings;

my %gene_name_to_gene_id;
{
    open(my $fh, "gene_ids.list") or die $!;
    while(<$fh>) {
        chomp;
        my $gene_id = $_;
        my ($gene_symbol, $rest) = split(/\|/, $gene_id);
        $gene_name_to_gene_id{$gene_symbol} = $gene_id;
    }
}


my %parafam_to_gene_list;
{
    open(my $fh, "genefam_list.txt") or die $!;
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $fam_id = $x[1];
        my $gene_symbol = $x[3];
        if (exists($gene_name_to_gene_id{$gene_symbol})) {
            push(@{$parafam_to_gene_list{$fam_id}}, $gene_symbol);
        }
    }
    close $fh;

}



foreach my $fam_id (keys %parafam_to_gene_list) {
    
    my @gene_symbols = @{$parafam_to_gene_list{$fam_id}};
    
    if (length(@gene_symbols) < 2) {
        next;
    }

    my $gene_id_1 = $gene_name_to_gene_id{$gene_symbols[0]};
    my $gene_id_2 = $gene_name_to_gene_id{$gene_symbols[1]};

    print "$fam_id\t$gene_id_1--$gene_id_2\n";
}

