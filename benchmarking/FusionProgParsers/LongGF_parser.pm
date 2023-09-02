package LongGF_parser;

use strict;
use warnings;
use Carp;

=LongGF_format

(lines start with SumGF and are space-delimited)

0 SumGF
1 ANAPC16:NFU1
2 1
3 chr10:72232999
4 chr2:69406021


=cut


sub parse_fusion_result_file {
    my ($file) = @_;

    open(my $fh, $file) or die "Error, cannot open file: $file";
    
    my @fusions;
    
    while (my $line = <$fh>) {
        
        if ($line !~ /^SumGF\s/) {
            next;
        }
        chomp $line;

        my @vals = split(/\s+/, $line);
        
        my $fusion_name = $vals[1];
        
        my ($fusion_gene_A, $fusion_gene_B) = split(/:/, $fusion_name);

        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions

        my $chr_coords_A = $vals[3];
        my $chr_coords_B = $vals[4];

        my ($chrA, $coordA) = split(/:/, $chr_coords_A);
        my ($chrB, $coordB) = split(/:/, $chr_coords_B);

        my $num_reads = $vals[2];
        
        
        my $struct = {
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

