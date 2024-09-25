#!/bin/bash

set -ex

rm -rf ./_*
rm -f ./preds.*

rm -f ./fusion_result_file_listing.dat ./max_F1_summary.tsv ./max_F1_summary.barplot.pdf ./max_F1_summary.linepoint.pdf ./breakpoint_* ./mean_PR_AUC_* ./pipe.log ./all_combined_results.*

