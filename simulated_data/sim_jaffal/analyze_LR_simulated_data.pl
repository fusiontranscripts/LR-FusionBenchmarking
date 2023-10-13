#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Cwd;
use File::Basename;
use lib ("$FindBin::Bin/../../PerlLib");
use Pipeliner;
use Process_cmd;
use DelimParser;
use Data::Dumper;

unless ($ENV{FUSION_ANNOTATOR}) {

    if (-d "$ENV{HOME}/GITHUB/CTAT_FUSIONS/FusionAnnotator") {
        $ENV{FUSION_ANNOTATOR} = "~/GITHUB/CTAT_FUSIONS/FusionAnnotator";
    }
    else {
        die "Error, must set env var FUSION_ANNOTATOR to point to base dir of\n"
            . "      git clone https://github.com/FusionAnnotator/FusionAnnotator.git\n"
            . "      (after having installed it)  ";
    }
}

unless ($ENV{TRINITY_HOME}) {
    die "Error, must specify env var TRINITY_HOME to trinity base installation directory";
}


my $usage = "\n\n\tusage: $0 [restrict_progs.file]\n\n";


my $benchmark_data_basedir = "$FindBin::Bin/../../";
my $benchmark_toolkit_basedir = "$FindBin::Bin/../../benchmarking";
my $fusion_annotator_basedir = $ENV{FUSION_ANNOTATOR};
my $trinity_home = $ENV{TRINITY_HOME};


my $sim_truth_set = "$FindBin::Bin/jaffal_sim_truth_set.tsv";

my $restrict_progs_file = $ARGV[0] || "";


main: {

    my $pipeliner = &init_pipeliner();
    
    ## create file listing
    my $cmd = "find ./prog_results -type f | ./util/make_LR_file_listing_input_table.pl $restrict_progs_file > fusion_result_file_listing.dat";
    $pipeliner->add_commands(new Command($cmd, "fusion_file_listing.ok"));
    
    # collect predictions
    $cmd = "./util/collect_LR_preds.pl fusion_result_file_listing.dat > preds.collected";
    $pipeliner->add_commands(new Command($cmd, "collect_preds.ok"));

    # map fusion predictions to gencode gene symbols based on identifiers or chromosomal coordinates.
    $cmd = "$benchmark_toolkit_basedir/map_gene_symbols_to_gencode.pl "
        . " preds.collected "
        . " $benchmark_data_basedir/resources/genes.coords.gz "
        . " > preds.collected.gencode_mapped ";

    $pipeliner->add_commands(new Command($cmd, "gencode_mapped.ok"));

    # annotate
    $cmd = "$fusion_annotator_basedir/FusionAnnotator --annotate preds.collected.gencode_mapped  -C 2 > preds.collected.gencode_mapped.wAnnot";
    $pipeliner->add_commands(new Command($cmd, "annotate_fusions.ok"));

    # filter HLA and mitochondrial features
    $cmd = "$benchmark_toolkit_basedir/filter_collected_preds.pl preds.collected.gencode_mapped.wAnnot > preds.collected.gencode_mapped.wAnnot.filt";
    $pipeliner->add_commands(new Command($cmd, "filter_fusion_annot.ok"));
    
    $pipeliner->run();
    
    ##################################
    ######  Scoring of fusions #######
    
    my $sample_to_truth_href = &parse_truth_set($sim_truth_set);
    
    my ($sample_to_fusion_preds_href, $column_headers_aref) = &parse_fusion_preds("preds.collected.gencode_mapped.wAnnot.filt");
    
    ########################################
    ## fusion scroring based on gene symbols, gene overlaps, and paralogs
    ## Breakpoints are not considered here.

    # score strictly
    &score_and_plot($sample_to_fusion_preds_href, 
                    $sample_to_truth_href, 
                    'analyze_strict', 
                    { allow_reverse_fusion => 0, allow_paralogs => 0 },
                    $column_headers_aref);
    
    # score allow reverse fusion
    &score_and_plot($sample_to_fusion_preds_href, 
                    $sample_to_truth_href, 
                    'analyze_allow_reverse', 
                    { allow_reverse_fusion => 1, allow_paralogs => 0 },
                    $column_headers_aref);

    # score allow reverse and allow for paralog-equivalence
    &score_and_plot($sample_to_fusion_preds_href, 
                    $sample_to_truth_href, 
                    'analyze_allow_rev_and_paralogs', 
                    { allow_reverse_fusion => 1, allow_paralogs => 1 },
                    $column_headers_aref);
    
    #############################
    ## fusion breakpoint analysis - only the breakpoint and distance to the breakpoint matters here.
    #############################
    
    &score_and_plot($sample_to_fusion_preds_href, 
                    $sample_to_truth_href, 
                    'breakpoint_exact', 
                    { breakpoint_eval => 1, # turns on breakpoint based evaluation.
                     max_dist => 0 },
                    $column_headers_aref);


    &score_and_plot($sample_to_fusion_preds_href, 
                    $sample_to_truth_href, 
                    'breakpoint_win10', 
                    { breakpoint_eval => 1,
                     max_dist => 10 },
                    $column_headers_aref);


    &score_and_plot($sample_to_fusion_preds_href, 
                    $sample_to_truth_href, 
                    'breakpoint_win100', 
                    { breakpoint_eval => 1,
                     max_dist => 100 },
                    $column_headers_aref);
    
    


    # generate summary max F1 plots:
    $cmd = "./util/plot_jaffal_F1_vs_divergence.Rscript";    
    $pipeliner->add_commands(new Command($cmd, "maxF1_summary_plots.ok"));

    # generate summary for breakpoint accuracy info:
    $cmd = "./util/plot_jaffal_breakpoint_accuracy.Rscript";
    $pipeliner->add_commands(new Command($cmd, "breakpoint_maxF1_summary_plots.ok"));

    $pipeliner->run();
    
    
    exit(0);
    
    
}


####
sub score_and_plot {
    my ($sample_to_fusion_preds_href, $truth_set_href, $analysis_token, $analysis_settings_href, $column_headers_aref) = @_;
    
    my $base_workdir = cwd();

    my $workdir = "__" . "$analysis_token";

    unless (-d $workdir) {
        mkdir ($workdir) or die "Error, cannot mkdir $workdir";
    }
    chdir ($workdir) or die "Error, cannot cd to $workdir";
    
    
    ####################################
    ## Examine each replicate separately
    
    #print Dumper($truth_set_href);
    #die;

    foreach my $sample_type (keys %$truth_set_href) {
        my $sample_checkpoint = "$sample_type.ok";
        if (! -e $sample_checkpoint) {
            &examine_sample($sample_type, $truth_set_href->{$sample_type}, $sample_to_fusion_preds_href->{$sample_type}, $analysis_settings_href, $column_headers_aref);
            &process_cmd("touch $sample_checkpoint");
        }
    }
    

    my $pipeliner = &init_pipeliner();

    my $cmd = 'find . -regex ".*fusion_preds.txt.scored.ROC"  >  ROC.files.list';
    $pipeliner->add_commands(new Command($cmd, "gather.ROC.files.ok"));

    $cmd = "../util/plot_jaffal_ROC_summary.Rscript";
    $pipeliner->add_commands(new Command($cmd, "plot_summary.ROC.files.ok"));


    $cmd = 'find . -regex ".*fusion_preds.txt.scored.PR.AUC"  >  PR.AUC.files.list';
    $pipeliner->add_commands(new Command($cmd, "gather.PR.AUC.files.ok"));
    
    
    $cmd = "../util/plot_jaffal_PR_AUC_barplot.Rscript";
    $pipeliner->add_commands(new Command($cmd, "plot_summary.PR.AUC.files.ok"));    
    


    $pipeliner->run();

    

    ######################################
    ## generate summary accuracy box plots
    
    # my $pipeliner = &init_pipeliner();
    # 
    # my $cmd = 'find . -regex ".*.scored.PR.AUC" -exec cat {} \\; > all.AUC.dat';
    # $pipeliner->add_commands(new Command($cmd, "gather_AUC.ok"));
    #
    # $cmd = "$benchmark_toolkit_basedir/plotters/AUC_boxplot.from_single_summary_AUC_file.Rscript all.AUC.dat";
    # $pipeliner->add_commands(new Command($cmd, "boxplot_rep_aucs.ok"));
    # 
    # $cmd = 'find . -regex ".*.scored" -exec cat {} \\; > all.scored.preds';
    # $pipeliner->add_commands(new Command($cmd, "gather_scores.ok"));
    # 
    # $pipeliner->run();
    
    
    # &ROC_and_PR("all.scored.preds");
        
    # examine sensitivity vs. expression level

    #if ($fusion_TPMs) {
    #    $cmd = "$benchmark_toolkit_basedir/fusion_preds_sensitivity_vs_expr.avg_replicates.pl all.scored.preds $fusion_TPMs > all.scored.preds.sensitivity_vs_expr.dat";
    #    $pipeliner->add_commands(new Command($cmd, "sens_vs_expr.avg_reps.ok"));
    #    
    #    $cmd = "$trinity_home/Analysis/DifferentialExpression/PtR  "
    #        . " -m all.scored.preds.sensitivity_vs_expr.dat "
    #        . " --heatmap "
    #        . " --sample_clust none --gene_clust ward "
    #        . " --heatmap_colorscheme 'black,purple,yellow'";
    #    $pipeliner->add_commands(new Command($cmd, "sens_expr_heatmap.ok"));
    #    
    #    $pipeliner->run();
    #}

    
    chdir $base_workdir or die "Error, cannot cd back to $base_workdir";
        
    return;
}

####
sub examine_sample {
    my ($sample_type, $sample_truth_href, $fusion_prog_to_preds_href, $analysis_settings_href, $column_headers_aref) = @_;
    
    my $basedir = cwd();

    my $sample_dir = "$sample_type";
    unless (-d $sample_dir) {
        mkdir($sample_dir) or die "Error, cannot mkdir $sample_dir";
    }
    chdir $sample_dir or die "Error, cannot cd to $sample_dir";

    my $sample_TP_fusions_file = "TP.fusions.list";
    my $fusion_preds_file = "fusion_preds.txt";

    my $prep_inputs_checkpoint = "_prep.ok";
    
    if (! -e $prep_inputs_checkpoint) {
        {
            my @TP_fusions = keys %$sample_truth_href;
            
            open (my $ofh, ">$sample_TP_fusions_file") or die "Error, cannot write to $sample_TP_fusions_file";
            print $ofh "fusion_name\tbreakpoint\tnum_reads\n";
            foreach my $fusion (@TP_fusions) {
                my $num_reads = $sample_truth_href->{$fusion}->{num_reads};
                my $breakpoint = $sample_truth_href->{$fusion}->{breakpoint};
                print $ofh join("\t", $fusion, $breakpoint, $num_reads) . "\n";
            }
            close $ofh;
        }
        
        {
            open (my $ofh, ">$fusion_preds_file") or die "Error, cannot write to $fusion_preds_file";
            my $delim_writer = new DelimParser::Writer($ofh, "\t", $column_headers_aref);
            
            foreach my $prog (keys %$fusion_prog_to_preds_href) {
                my @rows = @{$fusion_prog_to_preds_href->{$prog}};
                foreach my $row (@rows) {
                    $delim_writer->write_row($row);
                }
            }
            close $ofh;
        }
    
        &process_cmd("touch $prep_inputs_checkpoint");
    }
    
    ##################
    # score TP, FP, FN

    my $pipeliner = &init_pipeliner();
    
    my $cmd;

    if ($analysis_settings_href->{breakpoint_eval}) {
        
        my $max_dist = $analysis_settings_href->{max_dist};

        $cmd = "$benchmark_toolkit_basedir/fusion_breakpoints_to_TP_FP_FN.py --truth_fusions $sample_TP_fusions_file --pred_fusions $fusion_preds_file --max_dist $max_dist > $fusion_preds_file.scored";
        


    }
    else {
        
        $cmd = "$benchmark_toolkit_basedir/fusion_preds_to_TP_FP_FN.pl --truth_fusions $sample_TP_fusions_file --fusion_preds $fusion_preds_file";
        
        if ($analysis_settings_href->{allow_reverse_fusion}) {
            $cmd .= " --allow_reverse_fusion ";
        }
        if ($analysis_settings_href->{allow_paralogs}) {
            $cmd .= " --allow_paralogs $benchmark_data_basedir/resources/paralog_clusters.dat ";
        }
        
        $cmd .= " > $fusion_preds_file.scored";
        
        #print $cmd;
        #die;
    }
    
    $pipeliner->add_commands(new Command($cmd, "tp_fp_fn.ok"));
    
    $pipeliner->run();
    
        
    &ROC_and_PR("$fusion_preds_file.scored", $analysis_settings_href);
    
    
    
    chdir $basedir or die "Error, cannot cd back to $basedir";
        
    return;

}
 
   
####
sub ROC_and_PR {
    my ($preds_scored, $analysis_settings_href) = @_;

    ## run analysis pipeline
    my $pipeliner = &init_pipeliner();

    ##############
    # generate ROC
    
    my $ROC_script = "$benchmark_toolkit_basedir/all_TP_FP_FN_to_ROC.pl";

    if ($analysis_settings_href->{breakpoint_eval}) {
        $ROC_script = "$benchmark_toolkit_basedir/all_TP_FP_FN_to_ROC.for_brkpts.pl";
    }
    
    my $cmd = "$ROC_script $preds_scored > $preds_scored.ROC"; 
    $pipeliner->add_commands(new Command($cmd, "roc.ok"));
    
    # plot ROC
    $cmd = "$benchmark_toolkit_basedir/plotters/plot_ROC.Rscript $preds_scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "plot_roc.ok"));
        
    # plot F1
    $cmd = "$benchmark_toolkit_basedir/plotters/plot_F1_vs_min_frags.R $preds_scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "plot_F1_vs_min_frags.ok"));

    $cmd = "$benchmark_toolkit_basedir/plotters/plot_peak_F1_scatter.R $preds_scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "plot_peak_F1_scatter.ok"));
    
    # plot TP vs FP counts according to min frags per prog
    $cmd = "$benchmark_toolkit_basedir/plotters/plot_TP_FP_vs_minSum_per_prog.R $preds_scored.ROC";
    $pipeliner->add_commands(new Command($cmd, "sim_plot_TP_FP_vs_minFrags.ok"));
    
    
    ###################################
    # convert to Precision-Recall curve

    $cmd = "$benchmark_toolkit_basedir/calc_PR.py --in_ROC $preds_scored.ROC --out_PR $preds_scored.PR | sort -k2,2gr | tee $preds_scored.PR.AUC";
    $pipeliner->add_commands(new Command($cmd, "pr.ok"));

    # plot PR  curve
    # $cmd = "$benchmark_toolkit_basedir/plotters/plotPRcurves.R $preds_scored.PR $preds_scored.PR.plot.pdf";
    # $pipeliner->add_commands(new Command($cmd, "plot_pr.ok"));
    
    # plot AUC barplot
    $cmd = "$benchmark_toolkit_basedir/plotters/AUC_barplot.Rscript $preds_scored.PR.AUC";
    $pipeliner->add_commands(new Command($cmd, "plot_pr_auc_barplot.ok"));

    $pipeliner->run();

    #die "DEBUGGING";
    
    return;
        
}

sub parse_truth_set {
    my ($tp_fusions_file) = @_;

    my %sample_to_truth;

    open(my $fh, $tp_fusions_file) or die "Error, cannot open file $tp_fusions_file";
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    
    while(my $row = $delim_parser->get_row()) {
        my $sample_name = $row->{sample};
        my $fusion_name = $row->{FusionName};
        $row->{breakpoint} = $row->{Hg38_LeftBreakpoint} . "--" . $row->{Hg38_RightBreakpoint};
        $sample_to_truth{$sample_name}->{$fusion_name} = $row;
    }
    close $fh;

    return(\%sample_to_truth);

}


####
sub parse_fusion_preds {
    my ($preds_file, $preds_header_sref) = @_;

    open(my $fh, $preds_file) or die "Error, cannot open file: $preds_file";
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");

    my @column_headers = $delim_parser->get_column_headers();

    my %sample_to_preds;
    
    while(my $row = $delim_parser->get_row()) {
        
        my $sample_name = $row->{sample};
        my $prog = $row->{prog};
        push (@{$sample_to_preds{$sample_name}->{$prog}}, $row);
    }

    return(\%sample_to_preds, \@column_headers);
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

