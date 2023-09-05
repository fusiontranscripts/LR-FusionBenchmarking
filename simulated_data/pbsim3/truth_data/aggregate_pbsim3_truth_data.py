#!/usr/bin/env python3

import sys, os, re
import pandas as pd
import glob
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    full_df = None

    tsv_files = glob.glob("*.tsv")

    for tsv_file in tsv_files:
        dataset = tsv_file.split(".")[0]
        logging.info("Parsing {}".format(tsv_file))
        df = pd.read_csv(tsv_file, sep="\t")
        df.rename(
            columns={
                "LeftBreakPoint": "Hg38_LeftBreakpoint",
                "RightBreakPoint": "Hg38_RightBreakpoint",
                "sim_data_name": "sample",
                "num_LR": "num_reads",
            },
            inplace=True,
        )

        if full_df is None:
            full_df = df
        else:
            full_df = pd.concat([full_df, df])

    full_df = full_df[
        [
            "sample",
            "FusionName",
            "Hg38_LeftBreakpoint",
            "Hg38_RightBreakpoint",
            "num_reads",
        ]
    ]

    full_df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
