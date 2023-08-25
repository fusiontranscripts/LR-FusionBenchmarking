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




converter = get_lifter("hg19", "hg38")

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
        os.path.basename(fasta).replace(".fasta.gz", "")
        with gzip.open(fasta, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                desc = record.description.split()
                # example: LCORL|ENSG00000178177.10--CT47B1|ENSG00000236446.2
                one_gene_pair = re.search(r"(\S+)\|.*--(\S+)\|.*", desc[0])
                if one_gene_pair:
                    one_gene_pair = list(one_gene_pair.groups())
                else:
                    raise Exception
                # example: chr4:18022357-18022154,17963655-17963526[-]
                one_end = desc[2].split(":")
                # example: chrX:120007874-120007719,120006594-120006457[-]
                another_end = desc[3].split(":")

                one_end_coord = re.search(r"(?<=-)(\S+)\[", one_end[1].split(",")[-1])
                if one_end_coord:
                    one_end_coord = one_end_coord.groups()
                else:
                    raise Exception
                gold_standard_breakpoints_hg19.append(
                    [
                        one_end[0],
                        int(one_end_coord[0]),
                        another_end[0],
                        int(another_end[1].split(",")[0].split("-")[0]),
                        "--".join(one_gene_pair),
                    ]
                )
    # did not consider strand
    # do we consider strand ?
    gold_standard_breakpoints_hg19 = pd.DataFrame(gold_standard_breakpoints_hg19)
    print(gold_standard_breakpoints_hg19.head())
    gold_standard_breakpoints_hg38 = pd.concat(
        [
            gold_standard_breakpoints_hg19.apply(
                lambda x: ":".join(list(map(str, converter[x[0]][x[1]][0]))[:2]), axis=1
            ),
            gold_standard_breakpoints_hg19.apply(
                lambda x: ":".join(list(map(str, converter[x[2]][x[3]][0]))[:2]), axis=1
            ),
            gold_standard_breakpoints_hg19.iloc[:, -1],
        ],
        axis=1,
    )
    gold_standard_breakpoints_hg38.loc[
        :, "Breakpoint"
    ] = gold_standard_breakpoints_hg38.apply(lambda x: "--".join(x[:2]), axis=1)
    gold_standard_breakpoints_hg38.loc[
        :, "BreakpointLexSort"
    ] = gold_standard_breakpoints_hg38.apply(lambda x: "--".join(sorted(x[:2])), axis=1)
    gold_standard_breakpoints_hg38.loc[
        :, "FusionLexSort"
    ] = gold_standard_breakpoints_hg38.apply(
        lambda x: "--".join(sorted(x[4].upper().split("--"))), axis=1
    )
    print(gold_standard_breakpoints_hg38.loc[:, "FusionLexSort"].shape)
    print(gold_standard_breakpoints_hg38.loc[:, "FusionLexSort"].unique().shape)
    gold_standard_breakpoints_hg38 = gold_standard_breakpoints_hg38.set_axis(
        gold_standard_breakpoints_hg38.loc[:, "FusionLexSort"].str.upper(), axis="index"
    )
    print(gold_standard_breakpoints_hg38.head())
    gold_standard_breakpoints_hg38.columns = [
        "LeftBreakPoint",
        "RightBreakPoint",
        "FusionName",
        "Breakpoint",
        "BreakpointLexSort",
        "FusionLexSort",
    ]
    gold_standard_breakpoints_hg38.loc[
        :, "FusionName"
    ] = gold_standard_breakpoints_hg38.loc[:, "FusionName"].str.upper()
    return gold_standard_breakpoints_hg38

