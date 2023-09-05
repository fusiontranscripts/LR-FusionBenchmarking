package FusionSeeker_parser;

use strict;
use warnings;
use Carp;
use DelimParser;

=FusionSeeker_format

0       ID
1       Gene1
2       Gene2
3       NumSupp
4       Chrom1
5       Breakpoint1
6       Chrom2
7       Breakpoint2
8       SupportingReads

0       GF01
1       SLC10A1
2       NATD1
3       42
4       chr14
5       69786095
6       chr17
7       21244225
8       f0bcd4cf-5d1a-729d-c634-255e9306460d,0bc708dc-f2a5-436a-e98b-0b353eb9a548,11624683-4732-61df-d137-3a8acb2f71a4,17e51951-1545-9ed2-80bc-6368effcad41,1e7ffcca-e825-64db-c193-f74a9c79bf46,2487a021-4b0c-8ddf-3a21-d6ca7d0b40a3,250995e2-55e7-492d-2582-0e0a80fef44b,2a5f30fe-f294-908e-e40e-7debf260fe7e,2dd80384-4460-24a7-43da-be14d064dc09,33ea1a9c-46f4-c1b4-65f3-8fbd57a800fc,3af8d9c8-7bb6-3719-64d0-cb7807859b28,3c425fd3-a83e-8117-ec15-94ceb9fee5f3,3cbb2272-4004-0b0e-f855-e76767506b93,41448add-4fcd-59ab-1e8a-4b93a1ec374a,71483fc0-7253-2ff6-4795-7949d242a2ee,71dfee2c-c2ab-c650-cf13-f017a18d4259,82bf6392-da81-6486-747e-40668103e882,835ae101-fbcd-66ae-5549-89e4491ce4a1,a22dc8ef-535d-70bb-2a79-451443ee2f4d,a47eb1e2-5183-0c66-02d2-a784317e0975,ae659382-0fa7-af85-5dfe-28b544f544ee,b67ed79f-ea9c-6320-b41b-7a602b8fd66d,c80a7f74-d8f7-4c15-94dc-703813b3c11f,d580b9a7-b4ec-52b8-24e5-599264912663


=cut


sub parse_fusion_result_file {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file | ") or die "Error, cannot open file: $file";
    }
    else {
        open($fh, $file) or die "Error, cannot open file: $file";
    }
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");

    my @fusions;
    
    while(my $row = $delim_parser->get_row()) {
        
        my $fusion_gene_A = $row->{Gene1};
        my $fusion_gene_B = $row->{Gene2};
        
        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions


        my $chrA = $row->{Chrom1};
        my $coordA = $row->{Breakpoint1};
        
        my $chrB = $row->{Chrom2};
        my $coordB = $row->{Breakpoint2};
        
        my $num_reads = $row->{NumSupp};
                
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

