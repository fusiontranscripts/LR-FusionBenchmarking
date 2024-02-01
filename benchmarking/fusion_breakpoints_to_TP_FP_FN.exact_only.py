#!/usr/bin/env python3

import sys, os, re
import pandas as pd
import gzip
import intervaltree
from intervaltree import Interval
from typing import Any
from collections import defaultdict
import warnings
import numpy as np
import argparse



def main():

    parser = argparse.ArgumentParser(description="assign TP, FP, and FN pred_class for fusions based on breakpoints", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--truth_fusions", type=str, required=True, help="truth fusions")
    parser.add_argument("--pred_fusions", type=str, required=True, help="predicted fusions")

    args = parser.parse_args()

    truth_fusions = args.truth_fusions
    pred_fusions = args.pred_fusions


    truth_fusions_df = pd.read_csv(truth_fusions, sep="\t")
    truth_fusions_df['lexsort_breakpoint'] = truth_fusions_df['breakpoint'].apply(lambda x: "--".join(sorted(x.split("--"))))
    truth_fusions_df.rename(columns={'fusion_name' : 'truth_fusion_name',
                                     'breakpoint' : 'truth_breakpoint',
                                     'num_reads' : 'truth_num_reads'},
                            inplace=True)

    
    pred_fusions_df = pd.read_csv(pred_fusions, sep="\t")
    pred_fusions_df['lexsort_breakpoint'] = pred_fusions_df['breakpoint'].apply(lambda x: "--".join(sorted(x.split("--"))))


    ## should only be one sample type!
    assert len(pred_fusions_df['sample'].unique()) == 1, "Error, num samples != 1 "

    #must copy the truth set for each program to be analyzed separately so FNs show up in each case.
    all_truth_dfs = None
    progs = pred_fusions_df['prog'].unique()
    for prog in progs:
        prog_truth_df = truth_fusions_df.copy()
        prog_truth_df['prog'] = prog
                                       
        if all_truth_dfs is None:
            all_truth_dfs = prog_truth_df
        else:
            all_truth_dfs = pd.concat([all_truth_dfs, prog_truth_df])
    

    
    pred_fusions_df = pd.merge(all_truth_dfs, pred_fusions_df, on=['prog', 'lexsort_breakpoint'], how='outer')
    
    pred_fusions_df.sort_values(by=['prog', 'lexsort_breakpoint', 'num_reads'], ascending=[True, True, False], inplace=True)
    
    
    def assign_TP_FP_FN(df_slice):
        categories = list()
        for _, row in df_slice.iterrows():
            if pd.isnull(row['truth_fusion_name']):
                categories.append('FP')

            else:
                # truth fusion specified. either TP or FN depending on pred fusion status
                if pd.isnull(row['fusion']):
                    categories.append('FN')
                else:
                    categories.append('TP')
        
        # only score each breakpoint once.
        if len(categories) > 1:
            categories[1:] = ["NA_" + x for x in categories[1:] ]

        df_slice['pred_class'] = categories

        return df_slice



    pred_fusions_df = pred_fusions_df.groupby(['prog', 'lexsort_breakpoint']).apply(assign_TP_FP_FN)
    
  
    pred_fusions_df.to_csv(sys.stdout, sep="\t", index=False)


    sys.exit(0)


# methods based on Alvin's code:

def get_genes_to_breakpts(
        breakpoint_pairs, max_breakpoints_distance
    ) -> tuple[dict[Any, Any], dict[Any, Any]]:
        
    left_trees_dict: dict[Any, Any] = {}
    right_trees_dict: dict[Any, Any] = {}

    
    # chr12:52846197--chr17:9981884
    for index, breakpt in enumerate(breakpoint_pairs):
        if not re.search("^chr[^\\:]+:\\d+--chr[^\\:]+:\\d+$", breakpt):
            warnings.warn(f"{breakpt} lacks expected formatting. Skipping.", Warning)
            continue
            
        left, right = breakpt.split("--")
        #print(left); print(right)
        left = left.split(":")
        right = right.split(":")
        
        if left[0] not in left_trees_dict:
            left_trees_dict[left[0]] = intervaltree.IntervalTree()

        left_trees_dict[left[0]].add(
            Interval(
                int(left[1]) - max_breakpoints_distance // 2,
                int(left[1]) + max_breakpoints_distance // 2 + 1,
                (index, int(left[1])),
            )
        )

        if right[0] not in right_trees_dict:
            right_trees_dict[right[0]] = intervaltree.IntervalTree()

        right_trees_dict[right[0]].add(
            Interval(
                int(right[1]) - max_breakpoints_distance // 2,
                int(right[1]) + max_breakpoints_distance // 2 + 1,
                (index, int(right[1])),
            )
        )
            
    return left_trees_dict, right_trees_dict


def overlap_breakpoints(gold_standard_breakpts,
                        predicted_breakpts, 
                        max_breakpoints_distance=0):
    
    
    strictly_equal_breakpt = list(
        np.intersect1d(
            np.unique(gold_standard_breakpts),
            np.unique(predicted_breakpts),
        )
    )

    total_overlap_count = 0
    hit_paralogs_gold_pair = []
    recovered_not_found_predicted_pair = []

    not_found_gold_brkpts = np.setdiff1d(
        np.unique(gold_standard_breakpts),
        np.unique(predicted_breakpts),
    )
    not_found_prediction_brkpts = np.setdiff1d(
        np.unique(predicted_breakpts),
        np.unique(gold_standard_breakpts),
    )

  
    df = None
    
    if len(strictly_equal_breakpt) > 0:
        df = pd.DataFrame({ 'truth_brkpts' : strictly_equal_breakpt,
                            'pred_brkpts' : strictly_equal_breakpt,
                            'dist_left' : 0,
                            'dist_right' : 0,
                            'brkpt_match_type' : 'ExactMatched',
                          })
    if max_breakpoints_distance == 0:

        if len(not_found_gold_brkpts) > 0:
            not_found_gold_brkpts_df = pd.DataFrame({ 'truth_brkpts' : not_found_gold_brkpts,
                                                       'brkpt_match_type' : 'ExactUnmatched' })
            if df is None:
                df = not_found_gold_brkpts_df
            else:  
                df = pd.concat( [df, not_found_gold_brkpts_df])

        if len(not_found_prediction_brkpts) > 0:
            not_found_prediction_brkpts_df = pd.DataFrame({ 'pred_brkpts' : not_found_prediction_brkpts,
                                                            'brkpt_match_type' : 'ExactUnmatched' })
            if df is None:
                df = not_found_prediction_brkpts_df
            else:
                df = pd.concat([ df, not_found_prediction_brkpts_df ])

        return df

    # continue to search for additional entries within a max_breakpoints_distance allowance
    # build search tree for the unmatched predictions
    left_breakpt_tree, right_breakpt_tree = get_genes_to_breakpts(
        list(not_found_prediction_brkpts), max_breakpoints_distance
    )

    # search the unmatched gold brkpts using the unmatched prediction trees
    
    missing_gold_breakpts = []
    recovered_predicted_brkpts = []
    for gold_breakpt in not_found_gold_brkpts:
        left, right = gold_breakpt.split("--")
        left = left.split(":")
        right = right.split(":")

        if (right[0] not in right_breakpt_tree) or (
            left[0] not in left_breakpt_tree
        ):
            # NOTE: skip alternative contigs from liftover
            # If no chromosome matched, use infinity for maximum distances
            # then no fp_distance_list exists
            warnings.warn(
                f"Missing one or both chromosome pairs\
                in the prediction: {left[0]}, {right[0]}!",
                Warning,
            )
            missing_gold_breakpts.append(gold_breakpt)
            continue

        # left: [chr, start], gold standard left break point
        # left_overlapped_intervals: intervals from interval
        # trees, (start, end, data), data is the breakpoint pair index
        left_overlapped_intervals = list(
            left_breakpt_tree[left[0]].at(int(left[1]))
        )
        right_overlapped_intervals = list(
            right_breakpt_tree[right[0]].at(int(right[1]))
        )
        if (
            len(left_overlapped_intervals) > 0
            and len(right_overlapped_intervals) > 0
            and any(
                [
                    i.data[0] == j.data[0]
                    for i in left_overlapped_intervals
                    for j in right_overlapped_intervals
                ]
            )
        ):
            # Matched hit
            # if both left and right breakpoint overlaps with gold standard
            # and they belong to the same breakpoint pair
            local_intersected_breakpts = []
            local_intersected_intervals = []
            local_max_distances = []
            local_left_distances = []
            local_right_distances = []
            for i in left_overlapped_intervals:
                for j in right_overlapped_intervals:
                    if i.data[0] == j.data[0]:
                        local_intersected_breakpts.append(not_found_prediction_brkpts[ i.data[0] ])
                        local_intersected_intervals.append(
                            f"{left[0]}:{i.data[1]}--{right[0]}:{j.data[1]}"
                        )
                        total_overlap_count += 1
                        local_left_distances.append(abs(int(left[1]) - i.data[1]))
                        local_right_distances.append(abs(int(right[1]) - j.data[1]))

                        local_max_distances.append(
                            max(
                                [
                                    abs(int(left[1]) - i.data[1]),
                                    abs(int(right[1]) - j.data[1]),
                                ]
                            )
                        )
            if len(local_intersected_intervals) > 1:
                # NOTE: this condition never happen in the simulation
                warnings.warn(
                    "One gold standard pair match to multiple\
                    predicted breakpoints after window extension!" + str(local_intersected_intervals),
                    Warning,
                    stacklevel=2,
                )

            df_local = pd.DataFrame({ 'truth_brkpts' : gold_breakpt,
                                       'pred_brkpts' : local_intersected_breakpts,
                                       'left_distances' : local_left_distances,
                                        'right_distances' : local_right_distances,
                                        'max_distance' :  local_max_distances,
                                        'brkpt_match_type' : 'InexactMatched'} )
            
            recovered_predicted_brkpts.extend(local_intersected_breakpts)
            
            # append matched breakpoints to df
            if df is None:
                df = df_local
            else:
                df = pd.concat([df, df_local])

        else: # no overlaps found
            missing_gold_breakpts.append(gold_breakpt)

    missing_gold_brkpts_df = None
    if missing_gold_breakpts:
        missing_gold_breakpts_df = pd.DataFrame({ 'truth_brkpts' : missing_gold_breakpts,
                                                  'brkpt_match_type' : 'InexactUnmatched'} )
        
        if df is not None:
            df = pd.concat([df, missing_gold_breakpts_df])
        else:
            # really nothing matched
            df = missing_gold_breakpts_df

    remaining_unfound_prediction_breakpts = np.setdiff1d(
        np.unique(not_found_prediction_brkpts),
        np.unique(recovered_predicted_brkpts) )

    if len(remaining_unfound_prediction_breakpts) > 0:
        remaining_unfound_prediction_breakpts_df = pd.DataFrame({'pred_brkpts' : remaining_unfound_prediction_breakpts,
                                                                 'brkpt_match_type' : 'InexactUnmatched'} )
        if df is not None:
            df = pd.concat([df, remaining_unfound_prediction_breakpts_df])
        else:
            df = remaining_unfound_prediction_breakpts_df
            
    if total_overlap_count == 0:
        warnings.warn("No breakpoints found even with window extension!")

    return df




if __name__=='__main__':
    main()

