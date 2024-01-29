#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


my $usage = "\n\n\tusage: $0 preds.collected.byProg min_truth [extra_truth_preds_file]\n\n";

my $preds_collected = $ARGV[0] or die $usage;
my $min_truth = $ARGV[1] or die $usage;
my $extra_truth_preds_file = $ARGV[2];

main: {
    
    my %extra_true_preds;
    if ($extra_truth_preds_file) {
        &populate_extra_truth_preds(\%extra_true_preds, $extra_truth_preds_file);
    }
        
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
        
        my $proxy_fusion_name = $row->{proxy_fusion_name} or die "Error, no proxy fusion name reported";
        
        if ($prog_count >= $min_truth || $extra_true_preds{$proxy_fusion_name} ) {
            # write truth
            $truth_writer->write_row($row);
        }
        elsif ($prog_count > 1) {
            # write unsure
            $unsure_writer->write_row($row);
        }
        else {
            # write unique
            $unique_writer->write_row($row);
        }
    }
    
    close $fh;
    close $ofh_truth;
    close $ofh_unsure;
    close $ofh_unique;
    
    exit(0);
}



####
sub populate_extra_truth_preds {
    my ($extra_true_preds_href, $extra_truth_preds_file) = @_;

    open(my $fh, $extra_truth_preds_file) or die "Error, cannot open file: $extra_truth_preds_file";
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    while (my $row = $delim_parser->get_row()) {
        
        my $proxy_fusion_name = $row->{proxy_fusion_name} or die "Error, no proxy fusion name!";
        $extra_true_preds_href->{$proxy_fusion_name} = 1;
        
        # allow lex sorted version too
        my $lex_sorted_fusion_name = $row->{lex_ordered_fusion_name} or die "Error, missing lex_ordered_fusion_name for xtra ";
        $extra_true_preds_href->{$lex_sorted_fusion_name} = 1;
        
    }
    
    return;
}


