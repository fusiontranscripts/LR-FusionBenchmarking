package PBfusion_parser;

use strict;
use warnings;
use Carp;


=PBfusion_format
0       #chr1
1       start1
2       end1
3       chr2
4       start2
5       end2
6       id
7       score
8       strand1
9       strand2
10      info
11      extra

0       chr7
1       26206484
2       26206485
3       chr11
4       27806719
5       27806720
6       BP167
7       MEDIUM
8       +
9       +
10      RC=1;MD=0;GN=CBX3,CBX3P1;GI=ENSG00000122565.17,ENSG00000177447.6;GC=1,1;CL=FUSION:1;MQ=25,13;MI=0.903509,0.85034;SA=0,1;BP=1;EX=1,1;IN=0,0
11      RN=550da927-6b5b-65c4-90e3-16526c124559-ccs;ON=.;CB=.


=cut


sub parse_fusion_result_file {
    my ($file) = @_;

    open(my $fh, $file) or die "Error, cannot open file: $file";

    my @fusions;
    
    while(my $line = <$fh>) {
        
        if ($line =~ /^\#/) { next; }
        chomp $line;
        my @vals = split(/\t/, $line);
        
        my $chrA = $vals[0];
        my $orientA = $vals[8];
        my $coordA = ($orientA eq "+") ? $vals[2] : $vals[1];
        
        my $chrB = $vals[3];
        my $orientB = $vals[9];
        my $coordB = ($orientB eq "+") ? $vals[4] : $vals[5];
        
        my ($geneA, $geneB, $num_reads);
        my $info = $vals[10];
        if ($info =~ /RC=(\d+);.*;GN=([^;]+);/) {
            $num_reads = $1;
            my @genes = split(/,/, $2);
            if (scalar @genes > 2) {
                print STDERR "Error, not sure how to handle multiple gene pairs: @genes\n\n$line\n";
                next;
            }
            ($geneA, $geneB) = @genes;
        }
        
        my $struct = {
            geneA => $geneA,
            geneB => $geneB,
            
            breakpoint => "${chrA}:${coordA}--${chrB}:${coordB}",
            
            num_reads => $num_reads,
        };
        
        push (@fusions, $struct);

    }
    
    close $fh;

    return(@fusions);
}


1; #EOM

