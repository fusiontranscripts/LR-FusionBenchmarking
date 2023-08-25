#!/usr/bin/env python3


import os
from typing import Any
from collections import defaultdict
from dataclasses import dataclass
from dataclasses import field
import warnings
import numpy as np
import pandas as pd
import intervaltree  # type: ignore
from Bio import SeqIO  # type: ignore
from intervaltree import Interval  # type: ignore
from liftover import get_lifter  # type: ignore
import re
import glob
import sys
import gzip


def main():
    template_sequences = "/seq/RNASEQ/public_ftp/STAR_FUSION_PAPER/SupplementaryData/sim_reads/fusion_transcript_sequences/*.fasta.gz"

    gold_standard_brkpts = load_gold_standard(template_sequences)

    gold_standard_brkpts.to_csv("/seq/RNASEQ/public_ftp/STAR_FUSION_PAPER/SupplementaryData/sim_reads/fusion_transcript_sequences/fusion_breakpoint_summary.tsv", sep="\t", index=False)

    sys.exit(0)
    

"""
Example fasta header:
>LCORL|ENSG00000178177.10--CT47B1|ENSG00000236446.2 LCORL|ENST00000382224.1--CT47B1|ENST00000371311.3 chr4:18022357-18022154,17974508-17974443,17964672-17964593,17963655-17963526[-] chrX:120007874-120007719,120006594-120006457[-] FusedAt:480
"""


def load_gold_standard(template_sequences):
    """
    Load all the template fasta
    and return the gold standard break points

    Args
    --------
    :param template_sequences: `str`, fasta paths
    """
    gold_standard_breakpoints_hg19 = []
    for fasta in glob.glob(template_sequences):
        #print(fasta)
        dataset = os.path.basename(fasta).replace(".fasta.gz", "")
        with gzip.open(fasta, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                desc = record.description.split()
                # example: LCORL|ENSG00000178177.10--CT47B1|ENSG00000236446.2
                fusion_gene_pair = re.search(r"(\S+)\|.*--(\S+)\|.*", desc[0])
                if fusion_gene_pair:
                    fusion_gene_pair = list(fusion_gene_pair.groups())
                else:
                    raise Exception

                left_gene, right_gene = fusion_gene_pair

                ###########################
                ## get left breakpoint info
                # example: chr4:18022357-18022154,17963655-17963526[-]

                left_fusion_info = desc[2]
                left_breakpoint_chrom = left_fusion_info.split(":")[0]

                left_breakpoint_coord = re.search(r"-(\d+)\[", left_fusion_info.split(",")[-1])
                if left_breakpoint_coord:
                    left_breakpoint_coord = int(left_breakpoint_coord.groups()[0])
                else:
                    raise Exception

                ############################
                ## get right breakpoint info
                # example: chrX:120007874-120007719,120006594-120006457[-]

                right_breakpoint_info = desc[3]
                right_breakpoint_chrom = right_breakpoint_info.split(":")[0]
                right_breakpoint_coord = int(right_breakpoint_info.split(":")[1].split(",")[0].split("-")[0])

                
                
                gold_standard_breakpoints_hg19.append(
                    [
                        left_gene,
                        left_breakpoint_chrom,
                        left_breakpoint_coord,
                        right_gene,
                        right_breakpoint_chrom,
                        right_breakpoint_coord,
                        "--".join(fusion_gene_pair),
                        dataset,
                        record.description
                    ]
                )

    gold_standard_breakpoints_hg19 = pd.DataFrame(gold_standard_breakpoints_hg19)

    gold_standard_breakpoints_hg19.columns = ["LeftGene", "LeftChrom", "Hg19_LeftBreakpointCoord",
                                              "RightGene", "RightChrom", "Hg19_RightBreakpointCoord",
                                              "FusionName", "SourceDataset", "FastaDescription"]
                                             

    print(gold_standard_breakpoints_hg19.head())

    
    ## Convert to hg38 using liftover

    converter = get_lifter("hg19", "hg38")
    """
    >>> converter['chr1'][20304873]
    [('chr1', 19978380, '+')]
    """
        
    def run_hg19_to_hg38_liftover(chrom, coord):
        
        converted_coord_info = converter[chrom][coord][0]
        
        converted_chrom = converted_coord_info[0]
        converted_coord = converted_coord_info[1]

        converted_breakpoint = f"{converted_chrom}:{converted_coord}"
        
        return converted_breakpoint


    gold_standard_breakpoints_hg38 = gold_standard_breakpoints_hg19.copy()

    gold_standard_breakpoints_hg38['Hg38_LeftBreakpoint'] = gold_standard_breakpoints_hg38.apply(
        lambda x: run_hg19_to_hg38_liftover(x['LeftChrom'], x['Hg19_LeftBreakpointCoord']), axis=1)

    gold_standard_breakpoints_hg38['Hg38_RightBreakpoint'] = gold_standard_breakpoints_hg38.apply(
        lambda x: run_hg19_to_hg38_liftover(x['RightChrom'], x['Hg19_RightBreakpointCoord']), axis=1)
    
    print(gold_standard_breakpoints_hg38.head()) 
    
    return gold_standard_breakpoints_hg38




if __name__=='__main__':
    main()
