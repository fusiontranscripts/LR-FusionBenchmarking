#!/usr/bin/env python

import sys, os, re


coord_transfer = dict()

hg19_coords =  open("hg19.v19.ref_annot.gtf.gene_spans.coords.coords_that_transfer_ok").readlines()
hg38_coords = open("hg19.v19.ref_annot.gtf.gene_spans.coords.coords_that_transfer_ok.hg38_coords").readlines()

assert(len(hg19_coords) == len(hg38_coords))

hg19_coords = [x.rstrip() for x in hg19_coords]
hg38_coords = [x.rstrip() for x in hg38_coords]


for i in range(len(hg19_coords)):
    coord_transfer[ hg19_coords[i] ] = hg38_coords[i]


with open("hg19.v19.ref_annot.gtf.gene_spans") as fh:
    for line in fh:
        line = line.rstrip()
        (ensg_id, chrom, lend, rend, strand, gene_symbol, gene_type) = line.split("\t")
        coordset = chrom + ":" + lend + "-" + rend
        if coordset in coord_transfer:
            new_coordset = coord_transfer[coordset]
            new_chrom, coordpair = new_coordset.split(":")
            new_lend, new_rend = coordpair.split("-")

            # make coords format
            print("\t".join([gene_symbol, new_chrom, new_lend, new_rend, "hg19_gencode_v19_liftover", "0"]))


