#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Data::Dumper;
use Set::IntervalTree;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;

my $usage = "\n\n\tusage: $0 preds.collected genes.coords.gz [DEBUG]\n\n";

my $preds_collected_file = $ARGV[0] or die $usage;
my $genes_coords_file = $ARGV[1] or die $usage;
#my $genes_aliases_file = $ARGV[2] or die $usage;
my $DEBUG = $ARGV[2] || 0;

###############################################################################
### a few handy globals ;-)

my %GENE_ID_TO_GENE_STRUCTS;  # Gene_id returns list of all structs w/ that id

my $FEATURE_COUNTER;

my %CHR_COORDS_TO_GENE_STRUCTS;  # key on chr:lend-rend, return list of structs.


our %CHR_TO_ITREE;

my $EXCLUDE_LIST = "^(Y_RNA|snoU13|N/A)"; ## make regex, include | separated entries

my $MAX_GENE_ALIASES = 1000; ## should be generous enough, warding against weirdnesses
my $MAX_GENE_RANGE_SEARCH = 2e6;

my %PRIMARY_TARGET_ACCS;

my %GENE_ALIASES;


my %CACHED_GENE_MAPPINGS; # speed up by not repeating expensive work


##############################################################################



main: {

    print STDERR "-init interval trees\n" if $DEBUG;
    &init_interval_trees($genes_coords_file);

    #print STDERR "-init aliases\n" if $DEBUG;
    #&init_aliases($genes_aliases_file, \%GENE_ALIASES);

    
    open (my $fh, $preds_collected_file) or die "Error, cannot read file $preds_collected_file";
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");
    
    my @column_headers = $delim_parser->get_column_headers();


    my $delim_writer = new DelimParser::Writer(*STDOUT, "\t", [@column_headers, "mapped_gencode_A_gene_list", "mapped_gencode_B_gene_list"]);
    
    open (my $ofh_excluded, ">$preds_collected_file.excluded") or die $!;
    my $excluded_delim_writer = new DelimParser::Writer($ofh_excluded, "\t", [@column_headers, "reason_excluded"]);
    
    my %counts_of_excluded_entries;
    
    while (my $row = $delim_parser->get_row()) {
                
        my $fusion_name = $row->{fusion};
                
        my ($geneA, $geneB) = split(/--/, $fusion_name);

        
        if ($geneA =~ /$EXCLUDE_LIST/ || $geneB =~ /$EXCLUDE_LIST/) {
            # print STDERR "-excluding: " . join("\t", $row->{sample}, $row->{prog}, $fusion_name) . "\n";

            $row->{reason_excluded} = "missing representation in gencode for geneA or geneB.  geneA ($geneA) or  geneB ($geneB) is in exclude list: $EXCLUDE_LIST ";
            $counts_of_excluded_entries{ $row->{prog} } ++;
            
            next;
        }
        
        #print STDERR "\t** fusion: [$fusion]  \n";
        my $gencode_A_genes = &get_gencode_overlapping_genes($geneA);
        
        my $gencode_B_genes = &get_gencode_overlapping_genes($geneB);

        if ($gencode_A_genes eq "." || $gencode_B_genes eq ".") {
            #print STDERR "-skipping prediction due to excluding gene symbol pair: " . join("\t", $row->{sample}, $row->{prog}, $fusion_name) . "\n";
            $counts_of_excluded_entries{ $row->{prog} } ++;

            $row->{reason_excluded} = "missing representation in gencode for geneA or geneB.  geneA ($geneA) -> $gencode_A_genes and geneB ($geneB) -> $gencode_B_genes";
            $excluded_delim_writer->write_row($row);
            next;
        }
        
        
        $row->{mapped_gencode_A_gene_list} = $gencode_A_genes;
        $row->{mapped_gencode_B_gene_list} = $gencode_B_genes;
        
        $delim_writer->write_row($row);
                
    }
    close $fh;    


    print STDERR "Counts of excluded entries by prog: " . Dumper(\%counts_of_excluded_entries);

    print STDERR "-done gencode symbol mappings.\n";
    
    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;
    
    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd $cmd died with ret $ret";
    }

    return;
}


####
sub get_gencode_overlapping_genes {
    my ($gene_id) = @_;

    print "-get overlapping genes for: [$gene_id]\n" if $DEBUG;
    
    my @gencode_overlapping_genes = &find_overlapping_gencode_genes($gene_id);
    
    my $gencode_genes_text = (@gencode_overlapping_genes) ? join(",", sort @gencode_overlapping_genes) : ".";
    
    return($gencode_genes_text);
    
}

####
sub find_overlapping_gencode_genes {
    my ($gene_id_listing) = @_;

    $gene_id_listing =~ s/,\s*$//; # remove trailing comma as seen in some entries

    my %gencode_genes;

    foreach my $gene_id (split(/,/, $gene_id_listing) ) {

        my @mapped_genes;
        
        if (my $cached_gene_mappings_aref = $CACHED_GENE_MAPPINGS{$gene_id}) {
            @mapped_genes = @$cached_gene_mappings_aref;
            print "\t\t** using cached gene mappings for $gene_id as: @mapped_genes\n" if $DEBUG;
        }
        else {

            print "\t-examining [$gene_id] for overlaps\n" if $DEBUG;
            unless ($gene_id =~ /\w/) { next; }
            
            # retain original gene_id if recognized as primary
            if ($PRIMARY_TARGET_ACCS{$gene_id}) {
                $gencode_genes{$gene_id} = 1;
                print "\t\t** [$gene_id] is flagged as primary.\n" if $DEBUG;
            }
            
            # retain if we identify and recognize an alias for the id 
            my $alias = $GENE_ALIASES{$gene_id};
            if ($alias && $PRIMARY_TARGET_ACCS{$alias}) {
                $gencode_genes{$alias} = 1;
                print "\t\t** alias for [$gene_id] = [$alias] is found as primary.\n" if $DEBUG;
            }
            
            @mapped_genes = &__map_genes($gene_id);
            
            @mapped_genes = grep { $_ !~ /^ensg/i } @mapped_genes; #not reporting ENSG vals in mapping list.
            
            $CACHED_GENE_MAPPINGS{$gene_id} = [@mapped_genes];
        }
        
        if (scalar(@mapped_genes) <= $MAX_GENE_ALIASES) { 
            
            foreach my $mapped_gene_id (@mapped_genes) {
                
                $gencode_genes{$mapped_gene_id} = 1;
            }
        }
        else {
            print "\t\t\t !! too many mapped genes for [$gene_id]. Not using mappings here. (mappings included: @mapped_genes)\n" if $DEBUG;
        }
    }
    

    my @candidate_gencode_genes = keys %gencode_genes;

    return(@candidate_gencode_genes);
    
}


my %reported_missing_gene;

sub __map_genes {
    my ($gene_id) = @_;

    print STDERR "\t\t\t-mapping gene: $gene_id to overlaps\n" if $DEBUG;
    
    my $gene_structs_aref = $GENE_ID_TO_GENE_STRUCTS{$gene_id} || $GENE_ID_TO_GENE_STRUCTS{ lc $gene_id};

    unless (ref $gene_structs_aref) {
        unless ($reported_missing_gene{$gene_id}) {
            print STDERR "-warning, no gene stored for identifier: [$gene_id] - fusions including this symbol will be ignored\n";
            $reported_missing_gene{$gene_id} = 1;
        }
        return ();
    }
    

    my %overlapping_genes = ($gene_id => 1); # init with gene_id
    
    foreach my $gene_struct (@$gene_structs_aref) {

        my $chr = $gene_struct->{chr};
        my $lend = $gene_struct->{lend};
        my $rend = $gene_struct->{rend};
        
        print "\t\t\t[$gene_id] given coords: $chr:$lend-$rend ... examining overlaps\n" if $DEBUG;
        
        my $search_dist = $rend - $lend + 1;
        if ($search_dist > $MAX_GENE_RANGE_SEARCH) {
            print "\t\t\t\tsearch dist: $search_dist exceeds max allowed: $MAX_GENE_RANGE_SEARCH\n" if $DEBUG;
            next;
        }
        
        my $itree = $CHR_TO_ITREE{$chr} or die "Error, no itree for chr [$chr], " . Dumper($gene_struct);
        
        my $overlaps_aref = $itree->fetch($lend, $rend);
        
        #print STDERR "-overlapping features: " . Dumper($overlaps_aref);
        
        foreach my $feature_id_aref (@$overlaps_aref) {
            my $chr_span_token = $feature_id_aref->[0];

            my @all_gene_structs = @{$CHR_COORDS_TO_GENE_STRUCTS{"$chr_span_token"}};
            foreach my $overlap_gene_struct (@all_gene_structs) {

                my $overlap_gene_id = $overlap_gene_struct->{gene_id};
                                
                $overlapping_genes{$overlap_gene_id} = 1;
                
                print "\t\t\t\t$overlap_gene_id overlaps $chr_span_token\n" if $DEBUG;
            }
        }
        
    }
    
    return(keys %overlapping_genes);
}


####
sub add_gene_struct {
    my ($struct, $is_target) = @_;

    my $gene_id = $struct->{gene_id};
    my $chr = $struct->{chr};
    my $lend = $struct->{lend};
    my $rend = $struct->{rend};

    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend); # just to be safe

    if ($EXCLUDE_LIST) {
        if ($gene_id =~ /$EXCLUDE_LIST/) { return; }
    }

   
    $FEATURE_COUNTER++;

    $struct->{FEATURE_COUNTER} = $FEATURE_COUNTER;

    unless ($chr =~ /chr/) {
        $chr = "chr$chr";
    }

    my $itree = $CHR_TO_ITREE{$chr};
    unless ($itree) {
        $itree = $CHR_TO_ITREE{$chr} = Set::IntervalTree->new;
    }

    if ( ! exists $CHR_COORDS_TO_GENE_STRUCTS{"$chr:$lend-$rend"}) {
        $itree->insert(["$chr:$lend-$rend"], $lend, $rend);
    }
    
    push (@{$GENE_ID_TO_GENE_STRUCTS{$gene_id}}, $struct);

    push (@{$CHR_COORDS_TO_GENE_STRUCTS{"$chr:$lend-$rend"}}, $struct);
    
    #print STDERR "-adding " . Dumper($struct);

    return;
}

####
sub init_interval_trees {
    my ($gene_coords_file) = @_;
    
    print STDERR "-parsing annotation gene spans\n";

    unless ($gene_coords_file =~ /\.gz$/) {
        die "Error, gene_coords file $gene_coords_file is supposed to be gzipped";
    }
    
    open(my $fh, "gunzip -c $gene_coords_file | ") or die "Error, cannot open file: $gene_coords_file";
    my $header_line = <$fh>;
    unless ($header_line =~ /^gene_id/) {
        croak "Error, didn't parse header from $gene_coords_file";
    }
    while (<$fh>) {
        chomp;
        my ($gene_id, $chr, $lend, $rend, $file, $primary_target) = split(/\t/);

        if ($primary_target) {
            $PRIMARY_TARGET_ACCS{$gene_id} = 1;
        }
        
        
        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend); # just to be sure.

        my $struct = { chr => $chr,
                       lend => $lend,
                       rend => $rend,
                       gene_id => $gene_id,
                       file => $file,
        };
        
        &add_gene_struct($struct, $primary_target);
        
    }
    
    return;
}


####
sub overlaps {
    my ($struct_A, $struct_B) = @_;

    if ($struct_A->{lend} < $struct_B->{rend}
        &&
        $struct_A->{rend} > $struct_B->{lend}) {
        return(1);
    }
    else {
        return(0);
    }

}

####
sub init_aliases {
    my ($gene_aliases_file, $gene_aliases_href) = @_;

    open(my $fh, $gene_aliases_file) or die "Error, cannot open file $gene_aliases_file";
    while (<$fh>) {
        chomp;
        my ($gene_symbol, $gene_id) = split(/\t/);
        $gene_aliases_href->{$gene_id} = $gene_symbol;
    }
    close $fh;

    return;
}
