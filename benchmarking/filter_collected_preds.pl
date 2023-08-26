#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


my $usage = "\n\n\tusage: $0 preds.collected.wAnnot [min_reads=1]\n\n";

my $preds_file = $ARGV[0] or die $usage;
my $min_reads = $ARGV[1] || 1;

open (my $fh, $preds_file) or die "Error, cannot open file $preds_file";

my $delim_parser_reader = new DelimParser::Reader($fh, "\t");
my @column_headers = $delim_parser_reader->get_column_headers();

my $delim_writer = new DelimParser::Writer(*STDOUT, "\t", \@column_headers);

while(my $row = $delim_parser_reader->get_row()) {

    my $fusion_name = $row->{fusion};
    my $num_reads = $row->{num_reads};
    my $annot = $row->{annots};
    
    unless($num_reads >= $min_reads) {
        next;
    }
    

    if ($fusion_name =~ /(^HLA\-)|\-HLA\-/ 
        ||
        ($annot && 
         ($annot =~ /chrM:/i
          ||
          $annot =~ /NEIGHBOR/
          ||
          $annot =~ /BLAST/
          ||
          $annot =~ /GTEx|BodyMap|DGD_PARALOGS|HGNC_GENEFAM|Greger_Normal|Babiceanu_Normal|ConjoinG/
          ||
          $fusion_name =~ /IG[HKL].*--IG[HKL]/   
         )
        )
        ) 
    {
        next;
    }

    # passed
    $delim_writer->write_row($row);
}


exit(0);

