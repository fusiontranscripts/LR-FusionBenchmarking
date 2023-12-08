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


# retain cosmic in either orientation as long as reported in proper orientation at least once.

my %retain_cosmic_any_orient;

while(my $row = $delim_parser_reader->get_row()) {
    my $fusion_name = $row->{fusion};
    my $annot = $row->{annots}; 

    if ($annot =~ /cosmic/i) {
        my ($geneA, $geneB) = split(/--|::/, $fusion_name);
        $retain_cosmic_any_orient{$fusion_name} = 1;
        $retain_cosmic_any_orient{"$geneB--$geneA"} = 1;
        $retain_cosmic_any_orient{"$geneB\:\:$geneA"} = 1;
    }
}


open ($fh, $preds_file) or die "Error, cannot open file $preds_file";
$delim_parser_reader = new DelimParser::Reader($fh, "\t");
my $delim_writer = new DelimParser::Writer(*STDOUT, "\t", \@column_headers);

while(my $row = $delim_parser_reader->get_row()) {

    my $fusion_name = $row->{fusion};
    my $num_reads = $row->{num_reads};
    my $breakpoint = $row->{breakpoint};
    my $annot = $row->{annots};
    
    unless($num_reads >= $min_reads) {
        next;
    }
    
    if ( (! exists $retain_cosmic_any_orient{$fusion_name}) 
         
         && 

       (  
          $fusion_name =~ /(^HLA\-)|\-HLA\-/ 
          ||
        
          $breakpoint =~ /chrM:/
          
          ||
          
          (defined($annot) && 
           
           ($annot =~ /NEIGHBOR/
            ||
            $annot =~ /chrM\b/
            ||
            $annot =~ /BLAST/
            ||
            $annot =~ /GTEx|BodyMap|DGD_PARALOGS|HGNC_GENEFAM|Greger_Normal|Babiceanu_Normal|ConjoinG/
            ||
            $fusion_name =~ /IG[HKL].*--IG[HKL]/   
           )
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

