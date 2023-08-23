package FlairFusion_parser;

use strict;
use warnings;
use Carp;
use DelimParser;

=FlairFusion_format

0       #name
1       spanning reads
2       mapping score(1 is good)
3       seq agreement near breakpoint (1 is good)
4       avg frac of reads at loci in fusion
5       3' breakpoint
6       5' breakpoint

0       TENT5C--UACA
1       76
2       0.957
3       0.469
4       0.809
5       3'-UACA-chr15-70702281-0.784
6       5'-TENT5C-chr1-117622841-0.835


=cut


sub parse_fusion_result_file {
    my ($file) = @_;


    my $fh = open($file) or die "Error, cannot open file: $file";
    my $delim_parser = new DelimParser::Reader($fh, "\t");

    my @fusions;
    
    while(my $row = $delim_parser->get_row()) {

        my $fusion_name = $row->{"#name"};
                
        my ($fusion_gene_A, $fusion_gene_B) = split(/--/, $fusion_name);
        
        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions

        
        my $chr_coords_A = $row->{'5\' breakpoint'};
        my $chr_coords_B = $row->{'3\' breakpoint'};
        
        my @chr_coords_A_vals = split(/-/, $chr_coords_A);
        my $chrA = $chr_coords_A_vals[2];
        my $coordA = $chr_coords_A_vals[3];

        my @chr_coords_B_vals = split(/-/, $chr_coords_B);
        my $chrB = $chr_coords_B_vals[2];
        my $coordB = $chr_coords_B_vals[3];
        
        my $num_reads = $row->{spanning reads};
                
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

