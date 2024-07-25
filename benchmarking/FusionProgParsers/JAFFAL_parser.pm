package JAFFAL_parser;

use strict;
use warnings;
use Carp;
use DelimParser;

=JAFFAL_format

(note, actually comma-delimited)

0       sample
1       fusion genes
2       chrom1
3       base1
4       strand1
5       chrom2
6       base2
7       strand2
8       gap (kb)
9       spanning pairs
10      spanning reads
11      inframe
12      aligns
13      rearrangement
14      contig
15      contig break
16      classification
17      known

0       ONT_fus_sim_95err.fastq
1       TNK1:SF3B6
2       chr17
3       7383089
4       +
5       chr2
6       24067851
7       -
8       Inf
9       0
10      129
11      FALSE
12      TRUE
13      TRUE
14      040f0aea-ef98-0a56-4b92-0f81bb4e3885
15      289
16      HighConfidence
17      -


=cut


sub parse_fusion_result_file {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file | ") or die "Error, cannot open file $file";
    }
    else {
        open($fh, $file) or die "Error, cannot open file: $file";
    }



    my $delim_parser = new DelimParser::Reader($fh, ",");
    
    my @fusions;
    
    while(my $row = $delim_parser->get_row()) {
        
        my $fusion_name = $row->{'fusion genes'};
        my ($fusion_gene_A, $fusion_gene_B) = split(/:/, $fusion_name);
        
        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions

        my $chrA = $row->{chrom1};
        my $coordA = $row->{base1};

        my $chrB = $row->{chrom2};
        my $coordB = $row->{base2};
        
        my $num_reads = $row->{'spanning reads'};

        ## get sample encoding
        my $dataset = $row->{sample};

        
        #my $confidence = $row->{classification};
        #unless ($confidence eq "HighConfidence") { next; }
                
        
        if ($dataset =~ /Pac|ONT/ && $dataset =~ /\d+err/) {
            ## jaffal sim data customization here.
            
            my $seqtype;
            if ($dataset =~ /Pac/) {
                $seqtype = "Pac";
            }
            elsif ($dataset =~ /ONT/) {
                $seqtype = "ONT";
            }
            else {
                die "Error, cannot ascertain Pac or ONT from $dataset";
            }
            
            my $divergence_level;
            if ($dataset =~ /(\d+)err/) {
                $divergence_level = $1;
            }
            else {
                die "Error, cannot ascertain divergence level from $dataset";
            }
            
            $dataset = "${seqtype}_${divergence_level}";
        }
        
        
        my $struct = {
            
            dataset => $dataset,
            
            geneA => $fusion_gene_A,
            
            geneB => $fusion_gene_B,
            
            breakpoint => "${chrA}:${coordA}--${chrB}:${coordB}",
            
            num_reads => $num_reads,
        };
        
        push (@fusions, $struct);

    }
    
    close $fh;

    return(@fusions);
}


1; #EOM

