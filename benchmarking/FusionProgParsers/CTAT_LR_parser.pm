package CTAT_LR_parser;

use strict;
use warnings;
use Carp;
use DelimParser;

=CTAT_LR_fusion_format
0       #FusionName
1       num_LR
2       LeftGene
3       LeftLocalBreakpoint
4       LeftBreakpoint
5       RightGene
6       RightLocalBreakpoint
7       RightBreakpoint
8       SpliceType
9       LR_accessions
10      LR_FFPM
11      annots
12      max_LR_FFPM
13      frac_dom_iso
14      above_frac_dom_iso

0       PIK3CD-AS1--ZNF512B
1       67
2       PIK3CD-AS1
3       1389
4       chr1:9654198:-
5       ZNF512B
6       9562
7       chr20:63967523:-
8       ONLY_REF_SPLICE
9       766dc4e1-2248-8104-a17e-6b1b39433b46,91b5243f-8936-8c19-652f-0da2c00d4fd2,3fd66dae-a231-ebbf-04ca-4b2ffd175e72,06f3bcfd-aa0c-0986-079b-3abe74448e91,2cdc606e-1595-1f4e-c573-b0ac5c24a972,cf4fef11-6a88-9f34-b75a-5ebd6b36d242,3ba35b7e-c9a7-0041-d029-71c0eaf5b16c,f7f858ef-1c15-b290-398d-305271f94212,6616b3d6-8120-d259-7752-27c0811f39a4,c70d153e-c187-ad44-de83-0bc5ba83e732,223cfe2d-d67b-31da-5d54-64b30b7c15f3,9338becb-69f9-eee1-2dd1-648c64f42193,256e36c4-b3d0-c634-c09b-85bd0c7d5178,813e7fe3-86b3-08da-86ba-fcc69f344072,1c5777f0-27d1-da08-4720-05281a280006,cbc820bb-817b-8668-7b58-c9d52fed6b51,54768d2d-25d7-11df-178e-6103bc845579,eb6aa517-6a2d-86c7-1699-d679cc7e2046,14de402b-3083-e44f-13d6-f57358546789,0c536df1-6244-35e7-d258-82ef5097dedc,b62eafa4-e153-07d4-98cb-83d774776c65,526c2bce-ac08-c8c1-e0cc-4ee2711227b1,2b79a74f-d5d8-01dd-eb29-8f980cf4b0ab,053aac0f-2136-230e-6143-0638693ac3aa,f68cdda1-3bea-be65-5827-b4ae9969096f,f0b51bae-cb16-7adb-13fd-acea132786f4,608daedf-39f2-732e-37e2-53d2bd87b86e,10e7c85a-81a8-c0a7-8cb0-72f1c793fd5c,a79642f4-c187-cd1e-fcdd-d4c1bdb162d3,548c9ff4-9dc6-27cb-ec1d-42f202b7a528,17d197fb-c585-f70f-1ca3-06f823f7ad80,a0f0d453-f867-ff70-aa5e-19def7f44e84,b9e37d84-9b09-f1e6-ccc5-69b64af539d9,a05618b3-70e7-be99-31a8-1f80ec550471,d7328492-1896-3abe-1ef8-0f623490be99,d8b83037-2e97-ef9c-6842-49cdfa6af881,05da49de-60db-6380-6e1e-f33635d31aef,05c2ebee-6ffe-8c99-6fee-d5faee8deaea,bda618a8-fa1d-6dec-d45f-2d86c35d8268,82100c12-b577-d7ff-8e16-c5cfc3deaa11,ee7fbbcf-bb83-b60f-b9a5-acc7e7b195ee,4483a494-8603-e3a9-b11b-b89a6be07bd6,89b51271-d696-048e-2022-37aa1826d830,1dfb1050-19cf-fe82-1c49-15b9cd48b4a1,363fcbd4-d5de-840f-ee7a-87464038d2f6,4bbedb56-aef0-bbc4-d7b0-1f2a5e6761e8,1c40c44c-3ddf-549a-336b-2ab71b80f213,8df5f59e-a55a-d7c4-5387-27ae6e8318a7,08b0129b-945f-0441-72f7-59975a24ff30,57a128a1-babc-2c0c-dba5-511368d2eb5c,d42e633d-4a8f-3021-23d0-7128b61393cc,d9b94855-a233-3c1d-31ba-9fc3c6b79150,dbbdbc30-6aad-836f-ebfb-4f0a7386671e,e4814231-e496-1a6d-0e0e-f1b526248a90,511e6c79-c494-b85e-7d99-85c4d80d1c52,31dcb084-81c4-295f-a25c-c81d9dc36b04,cb15a66e-2d55-9e81-4c3a-7f78f6e92cc5,5c927976-e454-9fcc-5639-68ea6350c9ff,c540fc0f-35cc-4dd9-e463-898006638ccb,7469bdd6-99cb-9c60-d673-49383f5f16c6,f922a7b2-3825-b7d6-c66f-a63c4212a0a9,558106d8-90c6-1a1d-786a-2ec795364b45,ebe36cc2-fa57-5afc-080f-3577332cf43b,f37416a0-43bf-fe66-66c5-3ae70020662c,0745de41-bcf4-58b6-25ba-c498c5ed1f49,5ed04f2e-33fa-3b62-6778-aa6924769933,5a5c1b1f-d090-d093-1818-a04fe7b8015a
10      3745.318
11      ["INTERCHROMOSOMAL[chr1--chr20]"]
12      3745.318
13      1.0
14      True


=cut


sub parse_fusion_result_file {
    my ($file) = @_;


    my $fh = open($file) or die "Error, cannot open file: $file";
    my $delim_parser = new DelimParser::Reader($fh, "\t");

    my @fusions;
    
    while(my $row = $delim_parser->get_row()) {
        
        my $fusion_gene_A = $row->{LeftGene};
        my $fusion_gene_B = $row->{RightGene};
        
        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions

        my $splice_type = $row->{SpliceType};
        
        unless ($splice_type eq "ONLY_REF_SPLICE") { next; } # otherwise, too many assembly artifacts
        
        my $chr_coords_A = $row->{LeftBreakpoint};
        my $chr_coords_B = $row->{RightBreakpoint};
        
        my ($chrA, $coordA, $orientA) = split(/:/, $chr_coords_A);
        my ($chrB, $coordB, $orientB) = split(/:/, $chr_coords_B);

        my $num_reads = $row->{num_LR};
        
        
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

