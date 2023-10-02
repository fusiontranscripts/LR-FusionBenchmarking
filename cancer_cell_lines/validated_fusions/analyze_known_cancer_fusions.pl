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


my $truth_set = "$FindBin::Bin/validated_fusions.tsv";

my $restrict_progs_file = $ARGV[0] || "";


main: {

    my $pipeliner = &init_pipeliner();
    
    ##################################
    ######  Scoring of fusions #######
    
    my $sample_to_truth_href = &parse_truth_set($truth_set);

    my $cmd = "$benchmark_toolkit_basedir/filter_collected_preds.pl ../preds.collected.gencode_mapped.wAnnot 1 > preds.collected.gencode_mapped.wAnnot.filt.min1";
    &process_cmd($cmd);
    
    my ($sample_to_fusion_preds_href, $column_headers_aref) = &parse_fusion_preds("preds.collected.gencode_mapped.wAnnot.filt.min1");
    
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

    foreach my $sample_type (keys %$sample_to_fusion_preds_href) {
        my $sample_checkpoint = "$sample_type.ok";
        if (! -e $sample_checkpoint) {
            &examine_sample($sample_type, $truth_set_href->{$sample_type}, $sample_to_fusion_preds_href->{$sample_type}, $analysis_settings_href, $column_headers_aref);
            &process_cmd("touch $sample_checkpoint");
        }
    }

    my $cmd = "cat */fusion_preds.txt.scored | egrep '^(TP|FN)' > $analysis_token.TPs_n_FNs";
    &process_cmd($cmd);

    $cmd = "../scripts/eval_min_agree.pl $analysis_token.TPs_n_FNs 1 ../../progs_select.txt > $analysis_token.TPs.matrix";
    &process_cmd($cmd);
    
    
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
            print $ofh "fusion_name\tnum_reads\n";
            foreach my $fusion (@TP_fusions) {
                print $ofh join("\t", $fusion, "1") . "\n";
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
    
    my $cmd = "$benchmark_toolkit_basedir/fusion_preds_to_TP_FP_FN.pl --truth_fusions $sample_TP_fusions_file --fusion_preds $fusion_preds_file";
    
    if ($analysis_settings_href->{allow_reverse_fusion}) {
        $cmd .= " --allow_reverse_fusion ";
    }
    if ($analysis_settings_href->{allow_paralogs}) {
        $cmd .= " --allow_paralogs $benchmark_data_basedir/resources/paralog_clusters.dat ";
    }

    $cmd .= " > $fusion_preds_file.scored";


    $pipeliner->add_commands(new Command($cmd, "tp_fp_fn.ok"));
    
    $pipeliner->run();
    
    
    chdir $basedir or die "Error, cannot cd back to $basedir";

    return;

}
 

sub parse_truth_set {
    my ($tp_fusions_file) = @_;

    my %sample_to_truth;

    open(my $fh, $tp_fusions_file) or die "Error, cannot open file $tp_fusions_file";
    while(<$fh>) {
        unless (/\w/) { next; }
        chomp;
        my ($sample_name, $fusion_name) = split(/\|/);
        $sample_to_truth{$sample_name}->{$fusion_name} = 1;
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

