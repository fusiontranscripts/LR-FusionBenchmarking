#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


my $usage = "\n\n\tusage: $0 preds.collected.byProg min_truth\n\n";

my $preds_collected = $ARGV[0] or die $usage;
my $min_truth = $ARGV[1] or die $usage;

main: {
    
    open(my $fh, $preds_collected) or die "Error, cannot open file $preds_collected";

    my $delim_parser = new DelimParser::Reader($fh, "\t");
    my @column_names = $delim_parser->get_column_headers();


    my $out_basename = basename($preds_collected);
    
    open(my $ofh_truth, ">$out_basename.min_${min_truth}.truth_set") or die $!;
    my $truth_writer = new DelimParser::Writer($ofh_truth, "\t", \@column_names);

    open(my $ofh_unsure, ">$out_basename.min_${min_truth}.unsure_set") or die $!;
    my $unsure_writer = new DelimParser::Writer($ofh_unsure, "\t", \@column_names);
    
    open(my $ofh_unique, ">$out_basename.min_${min_truth}.unique_set") or die $!;
    my $unique_writer = new DelimParser::Writer($ofh_unique, "\t", \@column_names);
    
    while(my $row = $delim_parser->get_row()) {
        
        
        my $prog_count = $row->{num_progs};
        
        if ($prog_count >= $min_truth) {
            $truth_writer->write_row($row);
        }
        elsif ($prog_count > 1) {
            $unsure_writer->write_row($row);
        }
        else {
            $unique_writer->write_row($row);
        }
    }
    
    close $fh;
    close $ofh_truth;
    close $ofh_unsure;
    close $ofh_unique;
    
    exit(0);
}
