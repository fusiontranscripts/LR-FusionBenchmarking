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
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);



my $MIN_READ_SUPPORT = 1;

my $usage = <<__EOUSAGE__;


####################
#
# Optional:
#
#  --min_read_support <int>      minimum read support (default: $MIN_READ_SUPPORT)
#
#  --restrict_fuzzy_breakpoints    only consider fusions within +/- 5 bases of a reference exon boundary
#
#####################


__EOUSAGE__

    ;


my $help_flag;
my $RESTRICT_FUZZY_BREAKPOINTS = 0;



&GetOptions ( 'h' => \$help_flag,
              'min_read_support=i' => \$MIN_READ_SUPPORT,
              'restrict_fuzzy_breakpoints' => \$RESTRICT_FUZZY_BREAKPOINTS,
    );



if ($help_flag) {
    die $usage;
}


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


if (basename(cwd()) ne "SGNex_ONT") {
    die "Error, must run this while in the SGNex_ONT/ directory.";
}


my $benchmark_data_basedir = "$FindBin::Bin/..";
my $benchmark_toolkit_basedir = "$FindBin::Bin/../benchmarking";
my $fusion_annotator_basedir = $ENV{FUSION_ANNOTATOR};
my $trinity_home = $ENV{TRINITY_HOME};


main: {

    my $pipeliner = &init_pipeliner();
    
    ## create file listing
    my $cmd = "find ./prog_results -type f | ./util/make_LR_file_listing_input_table.pl > fusion_result_file_listing.dat";
    $pipeliner->add_commands(new Command($cmd, "fusion_file_listing.ok"));
    
    # collect predictions
    $cmd = "./util/collect_LR_preds.pl fusion_result_file_listing.dat > preds.collected";
    $pipeliner->add_commands(new Command($cmd, "collect_preds.ok"));


    # require reproducbily found fusions - at least two replicates predicted w/ fusion
    $cmd = "./util/extract_replicated_fusions.Rscript preds.collected > preds.collected.reproducible";
    $pipeliner->add_commands(new Command($cmd, "reproducible.ok"));

    my $preds_file_use = "preds.collected.reproducible";

    # map fusion predictions to gencode gene symbols based on identifiers or chromosomal coordinates.
    $cmd = "$benchmark_toolkit_basedir/map_gene_symbols_to_gencode.pl "
        . " $preds_file_use "
        . " $benchmark_data_basedir/resources/genes.coords.gz "
        . " > preds.collected.gencode_mapped ";
    
    $pipeliner->add_commands(new Command($cmd, "gencode_mapped.ok"));

    # annotate
    $cmd = "$fusion_annotator_basedir/FusionAnnotator --annotate preds.collected.gencode_mapped  -C 2 --include_reciprocal > preds.collected.gencode_mapped.wAnnot";
    $pipeliner->add_commands(new Command($cmd, "annotate_fusions.ok"));

    # filter HLA and mitochondrial features, and require min read support
    $cmd = "$benchmark_toolkit_basedir/filter_collected_preds.pl preds.collected.gencode_mapped.wAnnot $MIN_READ_SUPPORT > preds.collected.gencode_mapped.wAnnot.filt";
    $pipeliner->add_commands(new Command($cmd, "filter_fusion_annot.ok"));
    
    # filter out messy fusions (those containing genes predicted in fusions by multiple programs across multple samples
    $cmd = "$benchmark_toolkit_basedir/exclude_messy_fusions.pl  preds.collected.gencode_mapped.wAnnot.filt progs_select.txt 3 ";
    $pipeliner->add_commands(new Command($cmd, "filter_messy.ok"));
    
    # capture counts of progs agree: (also writes $preds_file.proxy_assignments ) with proxy fusion selected.
    my $cmd = "./util/SGNex_collected_preds_to_fusion_prog_support_listing.pl preds.collected.gencode_mapped.wAnnot.filt.pass progs_select.txt SGNEx-as_truth_fusions.lex_ordered.tsv > preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree";
    $pipeliner->add_commands(new Command($cmd, "byProgAgree.ok"));
        

    if ($RESTRICT_FUZZY_BREAKPOINTS) {
        # filter allow for fuzzy breakpoints only
        $cmd = "./util/filter_require_fuzzy_breakpoint.py preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments > preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.FUZZY_RESTRICTED";
        $pipeliner->add_commands(new Command($cmd, "fuzzy_filter.ok"));

        $preds_file_use = "$preds_file_use.fuzzy_ok";
    }
    
       

    
    $pipeliner->run();
    
    exit(0);
        
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

