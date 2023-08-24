#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use Carp;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

#################################################################################################
#
# Required:
#
#  --truth_fusions <string>   file containing a list of the true fusions.
#
#  --fusion_preds <string>    fusion predictions ranked accordingly.
#
#
# Optional:
#
#  --unsure_fusions <string>   fusions where we're not sure if it's a TP or FP
#
#  --allow_reverse_fusion     if true fusion is A--B, allow for it to be reported as B--A
#
#  --allow_paralogs <string>  file containing tab-delimited list of paralog clusters
#                             so if TP_A--TP_B  is a true fusion,
#                               paraA1--paraB2  would be considered an ok proxy and scored as a TP.
#
##################################################################################################


__EOUSAGE__


    ;



my $help_flag;
my $fusion_preds_file;
my $truth_fusions_file;

my $ALLOW_REVERSE_FUSION = 0;

my $ALLOW_PARALOGS = 0;
my $paralogs_file;
my $unsure_fusions_file;

&GetOptions ( 'h' => \$help_flag,
              'fusion_preds=s' => \$fusion_preds_file,
              'truth_fusions=s' => \$truth_fusions_file,
              'unsure_fusions=s' => \$unsure_fusions_file,
              
              'allow_reverse_fusion' => \$ALLOW_REVERSE_FUSION,
              
              'allow_paralogs=s' => \$paralogs_file,
    );


if ($help_flag) { die $usage; }

unless ($fusion_preds_file && $truth_fusions_file) {
    
    die $usage;

}

if (@ARGV) {
    die "Error, don't understand options: @ARGV";
}

if ($paralogs_file) {
    $ALLOW_PARALOGS = 1;
}

my %TP_fusions = &parse_fusion_listing($truth_fusions_file);

my %unsure_fusions;
if ($unsure_fusions_file) {
    %unsure_fusions = &parse_fusion_listing($unsure_fusions_file);
}

my %FP_progFusions;
my %seen_progTP;

my %paralog_fusion_to_TP_fusion;
if ($ALLOW_PARALOGS) {
    %paralog_fusion_to_TP_fusion = &parse_paralogs_integrate_parafusions(\%TP_fusions, $paralogs_file);
}



main : {

    my %prog_names;
    
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    my @column_headers = $delim_parser->get_column_headers();
    
    
    my $delim_writer = new DelimParser::Writer(*STDOUT, "\t", ['pred_class', @column_headers, 'selected_fusion', 'explanation']);
    
    my $sample_name;

    while (my $row = $delim_parser->get_row()) {
        
        if ($sample_name) {
            if ($sample_name ne $row->{sample}) {
                confess "Error, mixed sample names in predictions file: $sample_name vs. $row->{sample}";
            }
        }
        else {
            $sample_name = $row->{sample};
        }
        
        my $prog_name = $row->{prog};
        my $fusion_name = $row->{fusion};
        my $num_reads = $row->{num_reads};
        my $mapped_A_list = $row->{mapped_gencode_A_gene_list};
        my $mapped_B_list = $row->{mapped_gencode_B_gene_list};
        
                
        ## ensure everything is being compared in a case-insensitive manner.
        $fusion_name = uc $fusion_name;

        $mapped_A_list = uc $mapped_A_list;
        $mapped_B_list = uc $mapped_B_list;
        
        $prog_names{$prog_name} = 1;
        
        my ($geneA, $geneB) = split(/--/, $fusion_name);
        
        my @partnersA = ($geneA);
        my @partnersB = ($geneB);
        
        foreach my $ele (split(/,/, $mapped_A_list)) {
            if ($ele && ! grep { $_ eq $ele } @partnersA) {
                push (@partnersA, $ele);
            }
        }
        foreach my $ele (split(/,/, $mapped_B_list)) {
            if ($ele && ! grep { $_ eq $ele } @partnersB) {
                push (@partnersB, $ele);
            }
        }
        
        my ($pred_result, $explanation, $fusion_selected) = &classify_fusion_prediction($prog_name, \@partnersA, \@partnersB);

        unless ($fusion_selected) {
            $fusion_selected = ".";
        }
    
        $row->{pred_class} = $pred_result;
        $row->{selected_fusion} = $fusion_selected;
        $row->{explanation} = $explanation;

        $delim_writer->write_row($row);
                        
    }

    
    ## Report false-negatives (known fusions not predicted)
    
    foreach my $prog_name (keys %prog_names) {
        foreach my $fusion_name (keys %TP_fusions) {
            if (! $seen_progTP{"$prog_name,$fusion_name"}) {

                my ($geneA, $geneB) = &decode_fusion($fusion_name);
                my $core_fusion_name = join("--", $geneA, $geneB);
                
                my $row = { 'pred_class' => "FN",
                            'sample' => $sample_name,
                            'prog' => $prog_name,
                            'fusion' => $core_fusion_name,
                            'breakpoint' => "NA", # add this later.
                            'num_reads' => 0,
                            'mapped_gencode_A_gene_list' => '.',
                            'mapped_gencode_B_gene_list' => '.',
                            'annots' => '.',
                            'selected_fusion' => '.',
                            'explanation' => 'prediction_lacking' 
                };
                
                $delim_writer->write_row($row);
                
            }
        }
    }
    
    exit(0);
    
    
}


####
sub classify_fusion_prediction {
    my ($prog_name, $partnerA_aref, $partnerB_aref) = @_;
   
    my @fusion_candidates;

    my $primary_fusion_name;

    # build candidate fusion list
    foreach my $partnerA (@$partnerA_aref) {
        foreach my $partnerB (@$partnerB_aref) {
            
            my $fusion_candidate = &encode_fusion($partnerA, $partnerB);
            
            #note, the primary A--B will show up first in the list.
            unless ($primary_fusion_name) {
                $primary_fusion_name = $fusion_candidate;
            }
            
            push (@fusion_candidates, $fusion_candidate);

            if ($ALLOW_REVERSE_FUSION) {
                my $fusion_candidate = &encode_fusion($partnerB, $partnerA);
                push (@fusion_candidates, $fusion_candidate);
            }
            
        }
    }
    


    my ($accuracy_token, $accuracy_explanation, $fusion_selected); 

    foreach my $fusion_name (@fusion_candidates) {

        my $using_para_proxy = undef;
        
        if ($ALLOW_PARALOGS && exists $paralog_fusion_to_TP_fusion{$fusion_name}) {

            my $para_fusion_name = $paralog_fusion_to_TP_fusion{$fusion_name};

            if ($fusion_name ne $para_fusion_name) {
                $using_para_proxy = $fusion_name;
                # now set as fusion name to use in analysis below
                $fusion_name = $para_fusion_name;
            }
                        
        }
        
        my $prog_fusion = "$prog_name,$fusion_name";
        
        ############################
        ## Check for already seen TP
        
        if ($seen_progTP{$prog_fusion}) {
            $accuracy_token = "NA-TP";
            $accuracy_explanation = "already scored $prog_fusion as TP";
            if ($using_para_proxy) {
                $accuracy_explanation .= " (para of $using_para_proxy)";
            }
            last; 
        }

        ############################
        ## Check for already seen FP
        
        elsif ($FP_progFusions{$prog_fusion}) {
            $accuracy_token = "NA-FP";
            $accuracy_explanation = "already scored $prog_fusion as FP";
            if ($using_para_proxy) {
                $accuracy_explanation .= " (para of $using_para_proxy)";
            }     
            last;
        }
        
        ###############################
        ## Check for new TP recognition
        
        elsif ($TP_fusions{$fusion_name}) {
            $accuracy_token = "TP";
            $seen_progTP{$prog_fusion} = 1;
            $accuracy_explanation = "first encounter of TP $prog_fusion";
            $fusion_selected = $fusion_name;
            if ($using_para_proxy) {
                $accuracy_explanation .= " (para of $using_para_proxy)";
            }
            last;
        }

        ###########
        ## Check to see if we should ignore it
        elsif (%unsure_fusions && $unsure_fusions{$fusion_name}) {
            $accuracy_token = "NA-UNCLASS";
            $accuracy_explanation = "not classifying $fusion_name, in unsure list";
            if ($using_para_proxy) {
                $accuracy_explanation .= " (para of $using_para_proxy)";
            }
            ## Not short-circuiting here... might have a better classification for an alternative candidate arrangment.
        }
    }
    
    if ($accuracy_token) {
        if ($accuracy_explanation !~ /$primary_fusion_name/) {
            # include the primary name to faciliate further study of the comment.
            $accuracy_explanation .= " ($primary_fusion_name)";
        }
    }
     
    unless ($accuracy_token) {
        # must be a FP
        $accuracy_token = "FP";
        
        if (! defined $prog_name ) {
            confess "Error, no prog_name defined";
        }
        if (! defined $primary_fusion_name) {
            confess "Error, no primary fusion name defined";
        }
        
        my $prog_fusion = "$prog_name,$primary_fusion_name";
        $accuracy_explanation = "first encounter of FP fusion $prog_fusion";
        $FP_progFusions{$prog_fusion} = 1;
    }
    
    return($accuracy_token, $accuracy_explanation, $fusion_selected);
    
}



####
sub parse_fusion_listing {
    my ($fusions_file) = @_;

    my %fusions;
    
    open (my $fh, $fusions_file) or die "Error, cannot open file $fusions_file";
    
    my $delim_reader = new DelimParser::Reader($fh, "\t");
    
    while(my $row = $delim_reader->get_row()) {
        my $fusion = $row->{fusion} || $row->{fusion_name};
        $fusion = uc $fusion;
        $fusions{$fusion} = 1;
    }
    close $fh;
    
    return(%fusions);

}


####
sub parse_paralogs_integrate_parafusions {
    my ($orig_fusions_href, $paralogs_file) = @_;
    
    my %gene_to_para_list;
    {
        open (my $fh, $paralogs_file) or die $!;
        while (<$fh>) {
            chomp;
            
            my @x = split(/\s+/, uc $_);  ## case insensitive
            
            foreach my $gene (@x) {
                $gene_to_para_list{$gene} = \@x;
            }
        }
        close $fh;
    }

    
    my %paralog_fusion_to_orig_fusion;

    my @orig_fusions = keys %$orig_fusions_href;

    foreach my $orig_fusion_name (@orig_fusions) {
        my ($geneA, $geneB) = split(/--/, $orig_fusion_name);
        
        my @paraA = ($geneA);
        if (my $para_aref = $gene_to_para_list{$geneA}) {
            @paraA = @$para_aref;
        }
        my @paraB = ($geneB);
        if (my $para_aref = $gene_to_para_list{$geneB}) {
            @paraB = @$para_aref;
        }

        foreach my $gA (@paraA) {
            foreach my $gB (@paraB) {

                my $para_fusion = "$gA--$gB";
                
                $paralog_fusion_to_orig_fusion{$para_fusion} = $orig_fusion_name;
            }
        }
    }
    
    return(%paralog_fusion_to_orig_fusion);

}


####
sub encode_fusion {
    my ($geneA, $geneB) = @_;

    my $fusion_name = "$geneA--$geneB";

    return($fusion_name);
}

####
sub decode_fusion {
    my ($fusion_name) = @_;

    my ($geneA, $geneB) = split(/--/, $fusion_name);
    
    return($geneA, $geneB);

}

