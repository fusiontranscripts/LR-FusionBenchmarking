package PBfusion_parser_v4;

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


0       chr13
1       36174475
2       36174476
3       chr6
4       31811985
5       31811986
6       BP391
7       MEDIUM
8       -
9       -
10      RC=47;MD=0.394;GN=CCDC169,SOHLH2,HSPA1L;GI=ENSG00000242715.6,ENSG00000120669.14,ENSG00000204390.9;GC=45,47,47;CL=FUSION:47;MQ=60,60;MI=0.967678,0.940336;SA=45,2;BP=47;EX=47,47;IN=0,0;AC=chr13:36174475-/chr6:31811985-,chr13:36227346-/chr13:36202092-
11      RN=ff06541c-e40e-d59f-afa5-51aa0e93afa4-molecule,846225ed-2a7d-faa0-9903-799487741960-molecule,a0


=cut


sub parse_fusion_result_file {
    my ($file) = @_;

    # print STDERR "-pbfusion parsing file: $file\n";

    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file | ") or die "Error, cannot open file $file";
    }
    else {
        open($fh, $file) or die "Error, cannot open file: $file";
    }
    
    my %unique_fusions;

    my $num_fusions = 0;
    my $num_problem_fusions = 0;
    
    while(my $line = <$fh>) {
        
        if ($line =~ /^\#/) { next; }
        chomp $line;
        my @vals = split(/\t/, $line);
        
        $num_fusions += 1;

        #my $chrA = $vals[0];
        #my $orientA = $vals[8];
        #my $coordA = ($orientA eq "+") ? $vals[2] : $vals[1];
        
        #my $chrB = $vals[3];
        #my $orientB = $vals[9];
        #my $coordB = ($orientB eq "+") ? $vals[4] : $vals[5];
        
        #my ($geneA, $geneB, $num_reads);
        my $info = $vals[10];
        
        $info =~ /RC=(\d+);.*;GN=([^;]+);/;
        my $num_reads = $1;
        my @genes = split(/,/, $2);
        
        unless ($num_reads > 0 && scalar(@genes) > 1) {
            die "Error, not parsing multiple genes and read-supported fusions from: $info";
        }

        $info =~ /;AC=([^;]+)$/;
        my $breakpoint_info = $1 or die "Error, no breakpoint info from $info";
        my @breakpoints = split(m|,|, $breakpoint_info);
        
        unless (scalar(@breakpoints) == scalar(@genes)-1 ) {
            print STDERR "\n\nPBFv4 Error, diff number of breakpoints and genes parsed from: $info\n\n";
            $num_problem_fusions += 1;
            next;
        }
        

        my $breakpoint_parser = sub { 
            my ($breakpoint) = @_;
            if ($breakpoint =~ /^([^\:]+):(\d+)([\+\-])$/) {
                my $chrom = $1;
                my $pos = $2;
                $pos += 1;
                my $orient = $3;
                return($chrom, $pos, $orient);
            }
            else {
                die "Error, cannot parse breakpoint: $breakpoint";
            }
        };
        
        
        eval {
            
            
            for (my $i = 1; $i <= $#genes; $i++) {
    
                my $breakpoint = $breakpoints[$i-1];
                my ($breakpointA, $breakpointB) = split(m|/|, $breakpoint);
                
                my $geneA = $genes[$i-1];
                my ($chromA, $posA, $orientA) = &$breakpoint_parser($breakpointA);
                
                my $geneB = $genes[$i];
                my ($chromB, $posB, $orientB) = &$breakpoint_parser($breakpointB);
                
                my $struct = {
                    geneA => $geneA,
                    geneB => $geneB,
                    
                    breakpoint => "${chromA}:${posA}--${chromB}:${posB}",
                    
                    num_reads => $num_reads,
                };
            

                my $token = join("^^^", sort {$a cmp $b} ($geneA, $geneB, $breakpointA, $breakpointB));
                
                #print STDERR "TOKEN: $token\n";
                
                if ( (! exists $unique_fusions{$token} ) || ( $unique_fusions{$token} < $num_reads) ) {
                    
                    # in case the same exact fusion is reported multiple times, take the one with the greatest read support.
                    
                    $unique_fusions{$token} = $struct;
                }
                
            }
        };

        if ($@) {
            die "Error encountered $@ in parsing info $info";
        }
    }
    
    close $fh;


    my @fusions = values %unique_fusions;

    
    if ($num_fusions > 0) {
        my $frac_failed = sprintf("%.1f", $num_problem_fusions / $num_fusions * 100);
        print STDERR "\n\n** PBFUSION_v4 Parse stats: total fusion: $num_fusions, and $num_problem_fusions ignored due to multiple gene pairs = $frac_failed % for $file\n\n\n"; 
    }
    else {
        print STDERR "\n\n** WARNING: no fusions parsed for $file with pbfusion v4\n";
    }
    
    return(@fusions);
}


1; #EOM

