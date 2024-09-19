#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


my $usage = "\n\n\tusage: $0 fusion_preds.collected progs_to_consider.txt\n\n";

my $preds_file = $ARGV[0] or die $usage;
my $progs_to_consider_file = $ARGV[1] or die $usage;

main: {

    my %progs_to_consider;
    {
        open(my $fh, $progs_to_consider_file) or die "Error, cannot open file $progs_to_consider_file";
        while (<$fh>) {
            s/^\s+|\s+$//g;
            my $prog = $_;
            if ($prog =~ /^\#/) { next; }
            
            $progs_to_consider{$prog} = 1;
        }
        close $fh;
    }
    
    
    
    ## count fusion names.
    my %fusion_name_prog_counter;
    
    
    my @rows;
    open (my $fh, $preds_file) or die "Error, cannot open file $preds_file";
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    my @column_headers = $delim_parser->get_column_headers();

    while (my $row = $delim_parser->get_row()) {
        my $sample_name = $row->{sample};
        my $prog = $row->{prog};
        

        
        my $fusion_name = $row->{fusion};
        my $sample_fusion_name = "${sample_name}|${fusion_name}";
        $row->{sample_fusion_name} = $sample_fusion_name;
        
        push (@rows, $row);
    
        if ($progs_to_consider{$prog}) {
            $fusion_name_prog_counter{$sample_fusion_name}->{$prog} = 1;
        }
        
    }
    close $fh;
    
    ## define fusion proxy name 
    foreach my $row (@rows) {
        my $sample_name = $row->{sample};
        my $prog = $row->{prog};
        my $fusion_name = $row->{fusion};
        my $sample_fusion_name = $row->{sample_fusion_name};
        my ($geneA, $geneB) = split(/--/, $fusion_name);
        my $mapped_gencode_A_gene_list = $row->{mapped_gencode_A_gene_list};
        my $mapped_gencode_B_gene_list = $row->{mapped_gencode_B_gene_list};
        
        my $recip_sample_fusion_name = "$sample_name|$geneB--$geneA";
        
        my $count_progs_sample_fusion_name = &get_count_sample_fusion_name($sample_fusion_name, \%fusion_name_prog_counter);
        my $count_progs_recip_sample_fusion_name = &get_count_sample_fusion_name($recip_sample_fusion_name, \%fusion_name_prog_counter);

        my $proxy_fusion_name;
        my $proxy_fusion_type;

        
        # check tie
        if ($count_progs_sample_fusion_name == $count_progs_recip_sample_fusion_name) {
            $proxy_fusion_name = ($sample_fusion_name lt $recip_sample_fusion_name) ? $sample_fusion_name : $recip_sample_fusion_name;
            $proxy_fusion_type = "tie_lt";
        }
        
        # check incoming > recip
        elsif ($count_progs_sample_fusion_name > 1 && $count_progs_sample_fusion_name > $count_progs_recip_sample_fusion_name) {
            $proxy_fusion_name = $sample_fusion_name;
            $proxy_fusion_type = "dominant_choice";
        }
        # check recip dominant
        elsif ($count_progs_recip_sample_fusion_name > 1) {
            $proxy_fusion_name = $recip_sample_fusion_name;
            $proxy_fusion_type = "recip_selected";
        }

        # appears to be unique - check the overlapping genes for a recognizable fusion
        else {
            if (my $alt_fusion_name = &examine_overlapping_genes_for_fusion_name($sample_fusion_name, $sample_name, \%fusion_name_prog_counter, 
                                                                                 $mapped_gencode_A_gene_list, $mapped_gencode_B_gene_list) ) {
                $proxy_fusion_name = $alt_fusion_name;
                $proxy_fusion_type = "overlap_rep";
            }
            else {
                $proxy_fusion_name = $sample_fusion_name;
                $proxy_fusion_type = "orig_name";
            }
        }
        
        $row->{proxy_fusion_name} = $proxy_fusion_name;
        $row->{proxy_fusion_type} = $proxy_fusion_type;

    }
    
    # get prog to proxy fusion info:
    my $proxy_fusion_file = "$preds_file.proxy_assignments";
    open(my $ofh, ">$proxy_fusion_file") or die "Error, cannot write $proxy_fusion_file";
    my $tab_writer = new DelimParser::Writer($ofh, "\t", ["proxy_fusion_name", "proxy_fusion_type", @column_headers]);
    

    my %fusion_to_prog;
    foreach my $row (@rows) {
        my $proxy_fusion_name = $row->{proxy_fusion_name};
        my $prog = $row->{prog};

        if ($progs_to_consider{$prog}) {
            # count towards truth set definition
            $fusion_to_prog{$proxy_fusion_name}->{$prog} = 1;
        }
        
        $tab_writer->write_row($row);
    }
    close $ofh;
    
        
    my @fusion_structs;
    foreach my $fusion_name (reverse sort { scalar(keys %{$fusion_to_prog{$a}}) <=> scalar(keys %{$fusion_to_prog{$b}}) } keys %fusion_to_prog) {
        

        my @prognames = sort keys %{$fusion_to_prog{$fusion_name}};;
        
        my $num_progs = scalar(@prognames);
        
        push (@fusion_structs, { fusion_name => $fusion_name,
                                 prognames => \@prognames,
                                 count => $num_progs,
              } );

        
    }

    @fusion_structs = reverse sort {$a->{count} <=> $b->{count} } @fusion_structs;

    print join("\t", "proxy_fusion_name", "prog_names", "num_progs") . "\n";
    foreach my $fusion_struct (@fusion_structs) {
        
        print join("\t", $fusion_struct->{fusion_name}, 
                   join(",", @{$fusion_struct->{prognames}}),
                   $fusion_struct->{count}
                   ) . "\n";
                
    }
    
    exit(0);
}



####
sub get_count_sample_fusion_name {
    my ($sample_fusion_name, $fusion_name_prog_counter_href) = @_;

    if (exists $fusion_name_prog_counter_href->{$sample_fusion_name}) {
        my $num_progs = scalar(keys %{$fusion_name_prog_counter_href->{$sample_fusion_name}});
        return($num_progs);
    }
    else {
        return(0);
    }
}

####
sub examine_overlapping_genes_for_fusion_name {
    my ($sample_fusion_name, $sample_name, 
        $fusion_name_prog_counter_href,
        $mapped_gencode_A_gene_list,
        $mapped_gencode_B_gene_list) = @_;

    my $best_score = 0;
    my $best_name;

    my @A_symbols = split(/,/, $mapped_gencode_A_gene_list);
    my @B_symbols = split(/,/, $mapped_gencode_B_gene_list);

    foreach my $A_symbol (@A_symbols) {

        foreach my $B_symbol (@B_symbols) {
            
            my $sample_fusion_test_name = "$sample_name|$A_symbol--$B_symbol";
            my $recip_sample_fusion_test_name = "$sample_name|$B_symbol--$A_symbol";
            
            my $sample_fusion_test_name_count = &get_count_sample_fusion_name($sample_fusion_test_name, $fusion_name_prog_counter_href);

            if ($sample_fusion_test_name ne $sample_fusion_name && $sample_fusion_test_name_count > $best_score) {
                $best_score = $sample_fusion_test_name_count;
                $best_name = $sample_fusion_test_name;
            }

            my $recip_sample_fusion_test_name_count = &get_count_sample_fusion_name($recip_sample_fusion_test_name, $fusion_name_prog_counter_href);

            if ($recip_sample_fusion_test_name ne $sample_fusion_name && $recip_sample_fusion_test_name_count > $best_score) {
                $best_score = $recip_sample_fusion_test_name_count;
                $best_name = $recip_sample_fusion_test_name;
            }
        }
    }
    
    return($best_name);
}
