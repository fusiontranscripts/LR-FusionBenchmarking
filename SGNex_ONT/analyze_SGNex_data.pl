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



my $MIN_READ_SUPPORT = 1;

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


my $benchmark_data_basedir = "$FindBin::Bin/..";
my $benchmark_toolkit_basedir = "$FindBin::Bin/../benchmarking";
my $fusion_annotator_basedir = $ENV{FUSION_ANNOTATOR};
my $trinity_home = $ENV{TRINITY_HOME};


my $VALIDATED_FUSIONS_FILENAME = "validated_fusions.from_JAFFAL_supp_table3.txt";

main: {

    my $pipeliner = &init_pipeliner();
    
    ## create file listing
    my $cmd = "find ./prog_results -type f | ./util/make_LR_file_listing_input_table.pl  > fusion_result_file_listing.dat";
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
    $cmd = "$fusion_annotator_basedir/FusionAnnotator --annotate preds.collected.gencode_mapped  -C 2 --include_reciprocal > preds.collected.gencode_mapped.wAnnot";
    $pipeliner->add_commands(new Command($cmd, "annotate_fusions.ok"));

    # filter HLA and mitochondrial features, and require min read support
    $cmd = "$benchmark_toolkit_basedir/filter_collected_preds.pl preds.collected.gencode_mapped.wAnnot $MIN_READ_SUPPORT > preds.collected.gencode_mapped.wAnnot.filt";
    $pipeliner->add_commands(new Command($cmd, "filter_fusion_annot.ok"));
 
    # capture counts of progs agree: (also writes $preds_file.proxy_assignments ) with proxy fusion selected.
    $cmd = "./util/SGNex_collected_preds_to_fusion_prog_support_listing.pl  preds.collected.gencode_mapped.wAnnot.filt  progs_select.txt $VALIDATED_FUSIONS_FILENAME  > preds.collected.gencode_mapped.wAnnot.filt.proxy_assignments.byProgAgree";
    $pipeliner->add_commands(new Command($cmd, "byProgAgree.ok"));
    
    $pipeliner->run();
    

    exit(0);
    
=ladeda1
    
    ########################
    ## Evaluate predictions
    ########################
    
    
    my $input_file = "preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments";
    $input_file = &ensure_full_path($input_file);

    my $prog_agree_listing = "preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree";
    $prog_agree_listing = &ensure_full_path($prog_agree_listing);


    my ($sample_to_fusion_preds_href, $fusion_preds_column_headers_aref) = parse_fusion_preds($input_file);

    my ($sample_to_unique_fusion_preds_href, $unique_fusions_column_headers_aref) = parse_unique_fusion_preds($prog_agree_listing);
    

    # get list of the individual samples here
    my @samples;
    {

        my %s;
        
        open(my $fh, "fusion_result_file_listing.dat") or die $!;
        while(<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $sample_name = $x[1];
            $s{$sample_name} = 1;
        }
        @samples = keys %s;
    }
    
    
    
    
    #########################################
    ## examine each of the samples separately

    my %TRUTH_FUSIONS;
    {
        open(my $fh, $VALIDATED_FUSIONS_FILENAME) or die $!;
        while(<$fh>) {
            chomp;
            my $sample_fusion_name = $_;
            my ($sample, $fusion_name) = split(/\|/, $sample_fusion_name);
            push (@{$TRUTH_FUSIONS{$sample}}, $fusion_name);
        }
        close $fh;
    }

    
    

    my $base_workdir = cwd();


    
    foreach my $sample (@samples) {

        my $core_sample_name = (split(/_/, $sample))[1];
        
        my $truth_fusions_aref = $TRUTH_FUSIONS{$core_sample_name} or die "Error, no truth fusions defined for sample $sample , $core_sample_name";
        
        chdir $base_workdir or die $!;

        my $fusion_preds_aref = $sample_to_fusion_preds_href->{sample};
        
        &example_sample($sample, $truth_fusions_aref, $fusion_preds_aref, $fusion_preds_column_headers_aref);

    }

=cut

    exit(0);
    
}

=ladeda2

####
sub examine_sample {
    my ($sample, $truth_fusions_aref, $fusion_preds_aref, $fusion_pred_column_headers_aref, $unique_fusion_preds_aref) = @_;
    
    my $analysis_token = "eval-$sample";
    my $workdir = "__" . "$analysis_token";
    
    unless (-d $workdir) {
        mkdir ($workdir) or die "Error, cannot mkdir $workdir";
    }
    chdir ($workdir) or die "Error, cannot cd to $workdir";
    
    # creates two files:
    my $truth_set_fname = &ensure_full_path(basename($prog_agree_listing) . "validated.truth_set");
    my $unsure_set_fname = &ensure_full_path(basename($prog_agree_listing) . ".nonunique.unsure_set");


    my $sample_TP_fusions_file = "TP.fusions.list";
    my $fusion_preds_file = "fusion_preds.txt";

    my $prep_inputs_checkpoint = "_prep.ok";
    

    if (! -e $prep_inputs_checkpoint) {
        {
            my @TP_fusions = @$truth_fusions_aref;
            
            open (my $ofh, ">$truth_set_fname") or die "Error, cannot write to $truth_set_fname";
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
    
}

####
sub parse_unique_fusion_preds {
    my ($prog_agree_listing) = @_;

    my %sample_to_unique_fusion_preds;
    
    open(my $fh, $prog_agree_listing) or die $!;
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    my @column_headers = $delim_parser->get_column_headers();

    while(my $row = $delim_parser->get_row()) {
        my $sample_fusion_name = $row->{proxy_fusion_name};
        my $sample_name = split(/\|/, $sample_fusion_name)[0];
        my $num_progs = $row->{num_progs};
        if ($num_progs == 1) {
            push (@{$sample_to_preds{$sample_name}}, $row);
        }
    }

    return(\%sample_to_preds, \@column_headers);
}

    

####
sub parse_fusion_preds {
    my ($preds_file) = @_;

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

=cut
        
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

