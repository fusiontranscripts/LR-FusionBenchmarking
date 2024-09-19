#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Cwd;
use File::Basename;
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;
use Process_cmd;
use DelimParser;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;


my $benchmark_data_basedir = "$FindBin::Bin/..";
my $benchmark_toolkit_basedir = "$FindBin::Bin/../benchmarking";
my $fusion_annotator_basedir = $ENV{FUSION_ANNOTATOR};
my $trinity_home = $ENV{TRINITY_HOME};


my $SUPPORTED_TRUTH_SET = &ensure_full_path("SGNEx-as_truth_fusions.lex_ordered.tsv");

main: {
    
    
    my $input_file = "preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments";
    $input_file = &ensure_full_path($input_file);

    my $prog_agree_listing = "preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree";
    $prog_agree_listing = &ensure_full_path($prog_agree_listing);
        
    my $base_workdir = cwd();

    my $analysis_token = "valid_set_TP_uniq_FP";
    
    my $workdir = "__" . "$analysis_token";

    unless (-d $workdir) {
        mkdir ($workdir) or die "Error, cannot mkdir $workdir";
    }
    chdir ($workdir) or die "Error, cannot cd to $workdir";

    # creates two files:
    my $truth_set_fname = &ensure_full_path(basename($prog_agree_listing) . ".agree.truth_set");
    my $unsure_set_fname = &ensure_full_path(basename($prog_agree_listing) . ".nonunique.unsure_set");

    {
        my %truth_fusions;
        {
            open(my $fh, $SUPPORTED_TRUTH_SET) or die $!;
            my $delim_parser = new DelimParser::Reader($fh, "\t");
            while(my $row = $delim_parser->get_row()) {
                my $truth_fusion = $row->{'lex_sorted_fusion_name'} or die "Error, no lex_sorted_fusion_name found: " . Dumper($row);
                $truth_fusions{$truth_fusion} = 1;
            }
        }
        
        open(my $fh, $prog_agree_listing) or die $!;
        my $reader = new DelimParser::Reader($fh, "\t");
        my @column_headers = $reader->get_column_headers();

        open(my $truth_ofh, ">$truth_set_fname") or die $!;
        open(my $unsure_ofh, ">$unsure_set_fname") or die $!;
        my $truth_writer = new DelimParser::Writer($truth_ofh, "\t", \@column_headers);
        my $unsure_writer = new DelimParser::Writer($unsure_ofh, "\t", \@column_headers);

        while(my $row = $reader->get_row()) {
            my $proxy_fusion_name = $row->{proxy_fusion_name} or die "Error, no proxy_fusion_name specified: " . Dumper($row);
            my $num_progs = $row->{num_progs} or die "Error, num progs not specified";

            # also consider alt orientation
            my $alt_fusion = &make_alt_fusion_name($proxy_fusion_name);
            
            if ($truth_fusions{$proxy_fusion_name} || $truth_fusions{$alt_fusion}) {
                $truth_writer->write_row($row);
            }
            elsif ($num_progs > 1) {
                $unsure_writer->write_row($row);
            }
        }
    }
    
    my $pipeliner = &init_pipeliner();
    
    ## Examine accuracy by applying unsure and paralog-equiv options

    foreach my $settings_href ( #{ allow_paralogs => 0, unsure_fusions => undef },
                                #{ allow_paralogs => 1, unsure_fusions => undef },
                                { allow_paralogs => 0, unsure_fusions => $unsure_set_fname },
                                { allow_paralogs => 1, unsure_fusions => $unsure_set_fname } ) {

        &evaluate_predictions($input_file, $truth_set_fname, $settings_href);

    }


    print STDERR "Done.\n";
    exit(0);
    
}


####
sub make_alt_fusion_name {
    my ($sample_fusion_name) = @_;

    my ($sample, $fusion_name) = split(/\|/, $sample_fusion_name);

    my ($geneA, $geneB) = split(/--/, $fusion_name);

    return("$sample|$geneB--$geneA");
}


####
sub evaluate_predictions {
    my ($input_file, $min_agree_truth_set, $analysis_settings_href) = @_;

    my $output_filename = "eval_supported";
    my $checkpoint_token = "eval_supported";
    {
        my @analysis_token_pts;
        if ($analysis_settings_href->{allow_paralogs}) {
            push (@analysis_token_pts, "okPara");
        }
        if ($analysis_settings_href->{unsure_fusions}) {
            push (@analysis_token_pts, "ignoreUnsure");
        }
        if (@analysis_token_pts) {
            my $analysis_token = join("_", @analysis_token_pts);
            $output_filename .= ".$analysis_token";
            $checkpoint_token .= ".$analysis_token";
        }
    }
    $output_filename .= ".results";
    
    ## run analysis pipeline
    my $pipeliner = &init_pipeliner();

    ##################
    # score TP, FP, FN
    
    my $cmd = "$benchmark_toolkit_basedir/fusion_preds_to_TP_FP_FN.original.pl --truth_fusions $min_agree_truth_set --fusion_preds $input_file";

    $cmd .= " --allow_reverse_fusion "; # always do this here. Sim data shows it's important for some progs.
    
    if ($analysis_settings_href->{allow_paralogs}) {
        $cmd .= " --allow_paralogs $benchmark_data_basedir/resources/paralog_clusters.dat ";
    }
    
    if ($analysis_settings_href->{unsure_fusions}) {
        $cmd .= " --unsure_fusions " . $analysis_settings_href->{unsure_fusions};
    }

    $cmd .= " > $output_filename.scored";

    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.tp_fp_fn.ok"));

    ##############
    # generate ROC
    
    $cmd = "$benchmark_toolkit_basedir/all_TP_FP_FN_to_ROC.pl $output_filename.scored > $output_filename.scored.ROC"; 
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.roc.ok"));
    
    # plot ROC
    $cmd = "$benchmark_toolkit_basedir/plotters/plot_ROC.Rscript $output_filename.scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.plot_roc.ok"));


    # plot F1
    $cmd = "$benchmark_toolkit_basedir/plotters/plot_F1_vs_min_frags.R $output_filename.scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.plot_F1_vs_min_frags.ok"));

    $cmd = "$benchmark_toolkit_basedir/plotters/plot_peak_F1_scatter.R $output_filename.scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.plot_peak_F1_scatter.ok"));
    
    # plot TP vs FP counts according to min frags per prog
    $cmd = "$benchmark_toolkit_basedir/plotters/plot_TP_FP_vs_minSum_per_prog.R $output_filename.scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.plot_TP_FP_vs_minFrags.ok"));
    
                
    ###################################
    # convert to Precision-Recall curve
    
    $cmd = "$benchmark_toolkit_basedir/calc_PR.py --in_ROC $output_filename.scored.ROC --out_PR $output_filename.scored.PR | sort -k2,2gr | tee $output_filename.scored.PR.AUC";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.pr.ok"));
    
    # # plot PR curve
    $cmd = "$benchmark_toolkit_basedir/plotters/plotPRcurves.R $output_filename.scored.PR $output_filename.scored.PR.plot.pdf";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.plot_pr.ok"));
    
    # plot AUC barplot
    $cmd = "$benchmark_toolkit_basedir/plotters/AUC_barplot.Rscript $output_filename.scored.PR.AUC";
    $pipeliner->add_commands(new Command($cmd, "$checkpoint_token.plot_pr_auc_barplot.ok"));
    
    $pipeliner->run();

    return;
    
}


####
sub init_pipeliner {
    
    my $pipeliner = new Pipeliner(-verbose => 2, -cmds_log => 'pipe.log');
    my $checkpoint_dir = cwd() . "/_checkpoints";
    unless (-d $checkpoint_dir) {
        mkdir $checkpoint_dir or die "Error, cannot mkdir $checkpoint_dir";
    }
    $pipeliner->set_checkpoint_dir($checkpoint_dir);

    return($pipeliner);
}

