#!/usr/bin/env python3

import re
import glob
import gzip
import abc
import os
from typing import Any
from collections import defaultdict
from dataclasses import dataclass
from dataclasses import field
import warnings
import numpy as np
import pandas as pd
#import intervaltree  # type: ignore
#from Bio import SeqIO  # type: ignore
#from intervaltree import Interval  # type: ignore
#from liftover import get_lifter  # type: ignore


converter = get_lifter("hg19", "hg38")


# TODO: use pandera to check DataFrame col types
@dataclass
class BaseFusion(object):
    """Base class for benchmarking gene fusion

    :param output_result_file: input file from fusion software output path
    :type output_result_file: str

    """

    __metaclass__ = abc.ABCMeta

    output_result_file: str
    gold_standard: pd.DataFrame
    result: pd.DataFrame = pd.DataFrame()
    lower_threshold: float = 0
    upper_threshold: float = 100
    max_breakpoints_distance: int = (
        0  # expand predicted breakpoints 10bp around the center
    )
    sorted: bool = True
    fusion_sort_label: str = "LexSort"
    breakpt_sort_label: str = "LexSort"
    paralog_pairs_dict: dict[Any, Any] = field(default_factory=dict)
    paralog_pairs: list[Any] = field(default_factory=list)
    # NOTE: jaffal specific field
    jaffal_sample_id: str | None = None

    @property
    def reads_field(self):
        return "num_LR"

    def __post_init__(self):
        # TODO: add assert for breakpoints and fusionnames
        self.parse_results()
        self.clean_results()
        if self.sorted:
            self.fusion_sort_label = "LexSort"
            self.breakpt_sort_label = "LexSort"
        else:
            self.fusion_sort_label = "Name"
            self.breakpt_sort_label = ""

    @abc.abstractmethod
    def parse_results(self):
        """Method documentation"""

    @abc.abstractmethod
    def clean_results(self):
        """Method documentation"""

    @property
    def gold_standard_breakpts(self):
        return self._gold_standard_breakpts

    @gold_standard_breakpts.setter
    def gold_standard_breakpts(self, value):
        self._gold_standard_breakpts = value

    @property
    def gold_standard_fusionnames(self):
        return self._gold_standard_fusionnames

    @gold_standard_fusionnames.setter
    def gold_standard_fusionnames(self, value):
        self._gold_standard_fusionnames = value

    @property
    def predicted_breakpoints(self) -> np.ndarray:
        return self._predicted_breakpoints

    @predicted_breakpoints.setter
    def predicted_breakpoints(self, value) -> None:
        self._predicted_breakpoints = value

    @property
    def predicted_fusionnames(self):
        return self._predicted_fusionnames

    @predicted_fusionnames.setter
    def predicted_fusionnames(self, value):
        self._predicted_fusionnames = value

    @property
    def fusionname_metrics(self):
        overlap = self.overlap_fusionnames()

        if len(self.gold_standard_fusionnames) != 0:
            # NOTE: paralogs may match multiple gold standard
            sensitivity = len(np.unique(overlap)) / len(self.gold_standard_fusionnames)
        else:
            sensitivity = np.nan

        if len(self.predicted_fusionnames) != 0:
            # NOTE: paralogs may match multiple gold standard
            precision = len(np.unique(overlap)) / len(self.predicted_fusionnames)
        else:
            precision = np.nan

        if (
            np.isnan(sensitivity)
            or np.isnan(precision)
            or (sensitivity + precision == 0)
        ):
            F1 = np.nan
        else:
            F1 = 2 * sensitivity * precision / (sensitivity + precision)
        return (
            len(np.unique(overlap)),
            len(np.unique(self.predicted_fusionnames)) - len(np.unique(overlap)),
            sensitivity,
            precision,
            F1,
        )

    @property
    def breakpt_metrics(self):
        """Compute the breakpoints overlap metrics"""
        overlap = self.overlap_breakpoints()

        if len(self.gold_standard_breakpts) != 0:
            sensitivity = len(np.unique(overlap)) / len(self.gold_standard_breakpts)
        else:
            sensitivity = np.nan

        if len(self.predicted_breakpoints) != 0:
            precision = len(np.unique(overlap)) / len(self.predicted_breakpoints)
        else:
            precision = np.nan

        if (
            np.isnan(sensitivity)
            or np.isnan(precision)
            or (sensitivity + precision == 0)
        ):
            F1 = np.nan
        else:
            F1 = 2 * sensitivity * precision / (sensitivity + precision)
        return (
            len(np.unique(overlap)),
            len(np.unique(self.predicted_breakpoints)) - len(np.unique(overlap)),
            sensitivity,
            precision,
            F1,
        )

    def overlap_fusionnames(self, return_gold_standard_names=True):
        self.gold_standard_fusionnames = self.gold_standard[
            f"Fusion{self.fusion_sort_label}"
        ].unique()
        self.predicted_fusionnames = self.result[
            f"Fusion{self.fusion_sort_label}"
        ].unique()
        if self.paralog_pairs_dict == {}:
            return np.intersect1d(
                np.unique(self.gold_standard_fusionnames),
                np.unique(self.predicted_fusionnames),
            )

        # NOTE: match paralogs
        # CORRECTION: use the non-equal gold standard fusion
        not_found_pairs = np.setdiff1d(
            np.unique(self.gold_standard_fusionnames),
            np.unique(self.predicted_fusionnames),
        )

        not_found_predicted_pairs = np.setdiff1d(
            np.unique(self.predicted_fusionnames),
            np.unique(self.gold_standard_fusionnames),
        )

        paralog_hits = 0
        hit_paralogs_gold_pair = []
        recovered_not_found_predicted_pair = []

        # NOTE: for each gold standard that is not directly equal to prediction
        for not_found_pair in not_found_pairs:
            pair = not_found_pair.split("--")
            if len(pair) > 2:
                continue  # ignore chain genes like A--B--C--D...

            left, right = pair
            left_paralogs = [left]
            right_paralogs = [right]
            for a in self.paralog_pairs_dict[left]:
                left_paralogs += self.paralog_pairs[a]
            for b in self.paralog_pairs_dict[right]:
                right_paralogs += self.paralog_pairs[b]
            # there exist at least one pair of paralogs
            # overlapped with gold standard pairs
            if (
                len(
                    np.intersect1d(
                        [f"{i}--{j}" for i in left_paralogs for j in right_paralogs],
                        not_found_predicted_pairs,
                    )
                )
                > 0
            ):
                # pick a random first paralog pair for listing all true positive output
                hit_paralogs_gold_pair.append(not_found_pair)
                if (
                    len(
                        np.intersect1d(
                            [
                                f"{i}--{j}"
                                for i in left_paralogs
                                for j in right_paralogs
                            ],
                            not_found_predicted_pairs,
                        )
                    )
                    > 1
                ):
                    warnings.warn(
                        """multiple paralogs of one gold standard pair 
                         match to multiple predicted fusion pairs!""",
                        Warning,
                        stacklevel=2,
                    )
                    print(
                        "multiple",
                        np.intersect1d(
                            [
                                f"{i}--{j}"
                                for i in left_paralogs
                                for j in right_paralogs
                            ],
                            not_found_predicted_pairs,
                        ),
                    )
                    recovered_not_found_predicted_pair.append(
                        ";".join(
                            list(
                                np.intersect1d(
                                    [
                                        f"{i}--{j}"
                                        for i in left_paralogs
                                        for j in right_paralogs
                                    ],
                                    not_found_predicted_pairs,
                                )
                            )
                        )
                    )
                else:
                    recovered_not_found_predicted_pair.append(
                        np.intersect1d(
                            [
                                f"{i}--{j}"
                                for i in left_paralogs
                                for j in right_paralogs
                            ],
                            not_found_predicted_pairs,
                        )[0]
                    )
                paralog_hits += len(
                    np.intersect1d(
                        [f"{i}--{j}" for i in left_paralogs for j in right_paralogs],
                        not_found_predicted_pairs,
                    )
                )

        print(paralog_hits)
        if paralog_hits == 0:
            warnings.warn("No paralogs found!", Warning)

        if return_gold_standard_names:
            return (
                list(
                    np.intersect1d(
                        np.unique(self.gold_standard_fusionnames),
                        np.unique(self.predicted_fusionnames),
                    )
                )
                + hit_paralogs_gold_pair
            )
        return (
            list(
                np.intersect1d(
                    np.unique(self.gold_standard_fusionnames),
                    np.unique(self.predicted_fusionnames),
                )
            )
            + recovered_not_found_predicted_pair
        )

    def overlap_breakpoints(self, return_gold_standard_names: bool = True) -> list[str]:
        self.predicted_breakpoints = self.result[
            f"Breakpoint{self.breakpt_sort_label}"
        ].unique()
        self.gold_standard_breakpts = self.gold_standard[
            f"Breakpoint{self.breakpt_sort_label}"
        ].unique()

        strictly_equal_breakpt = list(
            np.intersect1d(
                np.unique(self.gold_standard_breakpts),
                np.unique(self.predicted_breakpoints),
            )
        )

        if self.max_breakpoints_distance == 0:
            return strictly_equal_breakpt

        total_overlap_count = 0
        hit_paralogs_gold_pair = []
        recovered_not_found_predicted_pair = []

        not_found_pair = np.setdiff1d(
            np.unique(self.gold_standard_breakpts),
            np.unique(self.predicted_breakpoints),
        )
        not_found_prediction_pair = np.setdiff1d(
            np.unique(self.predicted_breakpoints),
            np.unique(self.gold_standard_breakpts),
        )
        left_breakpt_tree, right_breakpt_tree = self.get_genes_to_breakpts(
            not_found_prediction_pair
        )

        for gold_breakpt in not_found_pair:
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
                local_intersected_intervals = []
                local_max_distances = []
                for i in left_overlapped_intervals:
                    for j in right_overlapped_intervals:
                        if i.data[0] == j.data[0]:
                            local_intersected_intervals.append(
                                f"{left[0]}:{i.data[1]}--{right[0]}:{j.data[1]}"
                            )
                            total_overlap_count += 1
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
                        predicted breakpoints after window extension!",
                        Warning,
                        stacklevel=2,
                    )
                    print(local_intersected_intervals)
                hit_paralogs_gold_pair.append(gold_breakpt)
                recovered_not_found_predicted_pair.append(
                    ";".join(local_intersected_intervals)
                )

                if recovered_not_found_predicted_pair[-1] == "":
                    print(
                        gold_breakpt,
                        left_overlapped_intervals,
                        right_overlapped_intervals,
                        local_intersected_intervals,
                    )
                    warnings.warn(
                        "Unmatched breakpoints for left and right end", Warning
                    )

        if total_overlap_count == 0:
            warnings.warn("No breakpoints found even with window extension!")

        if return_gold_standard_names:
            return strictly_equal_breakpt + hit_paralogs_gold_pair
        return strictly_equal_breakpt + recovered_not_found_predicted_pair

    def breakpoints_distances_to_goldstandard(self) -> tuple:
        """Return distances to the closest gold standards break points"""
        assert self.max_breakpoints_distance > 0
        not_found_pair = np.setdiff1d(
            np.unique(self.gold_standard_breakpts),
            np.unique(self.predicted_breakpoints),
        )
        not_found_prediction_pair = np.setdiff1d(
            np.unique(self.predicted_breakpoints),
            np.unique(self.gold_standard_breakpts),
        )
        left_breakpt_tree, right_breakpt_tree = self.get_genes_to_breakpts(
            not_found_prediction_pair
        )

        # only records distances when break points
        # are not strictly equal
        max_distances = []
        # these two variables contain the gold standard
        # lexsorted fusionnames or breakpoints keys
        # it does not have to be a match to the prediction in the window
        # used to measure all the closest breakpoints distance
        distance_keys = []
        breakpoints_keys = []

        for gold_breakpt in not_found_pair:
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
                local_intersected_intervals = []
                local_max_distances = []
                for i in left_overlapped_intervals:
                    for j in right_overlapped_intervals:
                        if i.data[0] == j.data[0]:
                            local_intersected_intervals.append(
                                f"{left[0]}:{i.data[1]}--{right[0]}:{j.data[1]}"
                            )
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
                        predicted breakpoints after window extension!",
                        Warning,
                        stacklevel=2,
                    )
                    print(local_intersected_intervals)

                # Maximum distances across multiple intersected intervals
                distance_keys.append(gold_breakpt)
                max_distances.append(max(local_max_distances))
                breakpoints_keys.append(";".join(local_intersected_intervals))
            else:
                # Measure the distance of closest breakpoint if not matched
                # even with window extensions on both ends
                max_distances_paired_breakpoints = np.array(
                    [
                        max(
                            [
                                abs(predicted_left_breakpt.data[1] - int(left[1])),
                                abs(predicted_right_breakpt.data[1] - int(right[1])),
                            ]
                        )
                        for predicted_left_breakpt in left_breakpt_tree[left[0]]
                        for predicted_right_breakpt in right_breakpt_tree[right[0]]
                        if predicted_left_breakpt.data[0]
                        == predicted_right_breakpt.data[0]
                    ]
                )
                paired_breakpoints = np.array(
                    [
                        f"{left[0]}:{predicted_left_breakpt.data[1]}--{right[0]}:{predicted_right_breakpt.data[1]}"
                        for predicted_left_breakpt in left_breakpt_tree[left[0]]
                        for predicted_right_breakpt in right_breakpt_tree[right[0]]
                        if predicted_left_breakpt.data[0]
                        == predicted_right_breakpt.data[0]
                    ]
                )
                # NOTE: no break pairs between two chromsome
                if len(paired_breakpoints) == 0:
                    warnings.warn(
                        f"No chromosome pairs in the \
                        prediction: {left[0]}, {right[0]}!",
                        Warning,
                    )
                else:
                    max_index = np.argmax(max_distances_paired_breakpoints)
                    max_distances.append(max_distances_paired_breakpoints[max_index])
                    distance_keys.append(gold_breakpt)
                    breakpoints_keys.append(paired_breakpoints[max_index])

        assert len(distance_keys) == len(breakpoints_keys)
        return distance_keys, breakpoints_keys, max_distances

    def get_genes_to_breakpts(
        self, breakpoint_pairs=None
    ) -> tuple[dict[Any, Any], dict[Any, Any]]:
        if breakpoint_pairs is None:
            breakpoint_pairs = self.result.BreakpointLexSort.values
        left_trees_dict: dict[Any, Any] = {}
        right_trees_dict: dict[Any, Any] = {}
        for index, breakpt in enumerate(breakpoint_pairs):
            left, right = breakpt.split("--")
            left = left.split(":")
            right = right.split(":")
            if left[0] not in left_trees_dict:
                left_trees_dict[left[0]] = intervaltree.IntervalTree()
                left_trees_dict[left[0]].add(
                    Interval(
                        int(left[1]) - self.max_breakpoints_distance // 2,
                        int(left[1]) + self.max_breakpoints_distance // 2,
                        (index, int(left[1])),
                    )
                )
            else:
                left_trees_dict[left[0]].add(
                    Interval(
                        int(left[1]) - self.max_breakpoints_distance // 2,
                        int(left[1]) + self.max_breakpoints_distance // 2,
                        (index, int(left[1])),
                    )
                )
            if right[0] not in right_trees_dict:
                right_trees_dict[right[0]] = intervaltree.IntervalTree()
                right_trees_dict[right[0]].add(
                    Interval(
                        int(right[1]) - self.max_breakpoints_distance // 2,
                        int(right[1]) + self.max_breakpoints_distance // 2,
                        (index, int(right[1])),
                    )
                )
            else:
                right_trees_dict[right[0]].add(
                    Interval(
                        int(right[1]) - self.max_breakpoints_distance // 2,
                        int(right[1]) + self.max_breakpoints_distance // 2,
                        (index, int(right[1])),
                    )
                )
        return left_trees_dict, right_trees_dict

    def analysis_FN_FP(self):
        gold_standard_template = self.gold_standard.copy()
        assert (
            len(self.gold_standard_breakpts)
            == np.unique(self.gold_standard_breakpts).shape[0]
        )

        # Search equal fusion names
        gold_standard_template.loc[
            :, f"Fusion{self.fusion_sort_label}Equal"
        ] = gold_standard_template.loc[:, f"Fusion{self.fusion_sort_label}"].isin(
            self.predicted_fusionnames
        )
        gold_standard_template.loc[
            :, f"Breakpoint{self.breakpt_sort_label}Equal"
        ] = gold_standard_template.loc[:, f"Breakpoint{self.breakpt_sort_label}"].isin(
            self.predicted_breakpoints
        )
        assert (
            gold_standard_template.index.shape[0]
            == gold_standard_template.index.unique().shape[0]
        ), "non unique index of FusionName"

        gold_standard_include_paralogs = self.overlap_fusionnames()
        prediction_include_paralogs = self.overlap_fusionnames(
            return_gold_standard_names=False
        )
        if (
            len(gold_standard_include_paralogs)
            != np.unique(gold_standard_include_paralogs).shape[0]
        ):
            warnings.warn(
                "Paralogs match of predicted fusion name \
                to multiple gold standard fusion name",
                Warning,
                stacklevel=2,
            )
            assert (
                len(gold_standard_include_paralogs)
                == np.unique(gold_standard_include_paralogs).shape[0]
            ), "multiple matching of paralogs to one gold standard found!"
        assert (
            len(prediction_include_paralogs)
            == np.unique(prediction_include_paralogs).shape[0]
        )
        gold_standard_template.loc[
            :, f"Fusion{self.fusion_sort_label}ParalogMatch"
        ] = False
        gold_standard_template.loc[
            gold_standard_include_paralogs,
            f"Fusion{self.fusion_sort_label}ParalogMatch",
        ] = True
        gold_standard_template.loc[
            gold_standard_include_paralogs,
            f"Fusion{self.fusion_sort_label}ParalogNames",
        ] = prediction_include_paralogs

        flatten_fusion_names = []
        for prediction in prediction_include_paralogs:
            if ";" in prediction:
                flatten_fusion_names += prediction.split(";")
            else:
                flatten_fusion_names.append(prediction)
        # False positive pairs that cannot directly match name or paralogs' names
        FP_genes = np.setdiff1d(self.predicted_fusionnames, flatten_fusion_names)

        gold_standard_include_extension = self.overlap_breakpoints()
        prediction_include_extension = self.overlap_breakpoints(
            return_gold_standard_names=False
        )
        (
            distance_keys,
            breakpoints_keys,
            max_distances,
        ) = self.breakpoints_distances_to_goldstandard()
        if (
            len(gold_standard_include_extension)
            != np.unique(gold_standard_include_extension).shape[0]
        ):
            warnings.warn(
                "Paralogs match of predicted fusion name \
                to multiple gold standard fusion name",
                Warning,
                stacklevel=2,
            )
            assert (
                len(gold_standard_include_extension)
                == np.unique(gold_standard_include_extension).shape[0]
            ), "multiple matching of paralogs to one gold standard found!"

        assert (
            len(self.predicted_breakpoints)
            == np.unique(self.predicted_breakpoints).shape[0]
        )
        assert (
            len(gold_standard_include_extension)
            == np.unique(gold_standard_include_extension).shape[0]
        )
        assert (
            len(prediction_include_extension)
            == np.unique(prediction_include_extension).shape[0]
        ), f"{len(prediction_include_extension)}\
        {np.unique(prediction_include_extension).shape}"
        assert len(distance_keys) == len(max_distances)
        # NOTE: after window extension, there are multiple
        # breakpoints matching one gold standard breakpoints
        print(
            len(prediction_include_extension),
            np.unique(prediction_include_extension).shape[0],
        )
        gold_standard_template = gold_standard_template.set_index(
            f"Breakpoint{self.breakpt_sort_label}"
        )
        gold_standard_template.loc[
            :, f"Breakpoint{self.breakpt_sort_label}WindowMatch"
        ] = False
        gold_standard_template.loc[
            :, f"Breakpoint{self.breakpt_sort_label}WindowNames"
        ] = ""
        gold_standard_template.loc[
            gold_standard_include_extension,
            f"Breakpoint{self.breakpt_sort_label}WindowMatch",
        ] = True
        gold_standard_template.loc[
            gold_standard_include_extension,
            f"Breakpoint{self.breakpt_sort_label}WindowNames",
        ] = prediction_include_extension
        # MaxDistance including every breakpoints that share chromosome pairs
        # It can be a overlapped breakpoints pairs maximum distance
        # Or pairs that not in the window
        gold_standard_template.loc[
            distance_keys, f"Breakpoint{self.breakpt_sort_label}MaxDistance"
        ] = max_distances
        gold_standard_template.loc[
            distance_keys,
            f"Breakpoint{self.breakpt_sort_label}MaxDistanceClosestBreakPoint",
        ] = breakpoints_keys
        gold_standard_template = gold_standard_template.set_index(
            f"Fusion{self.fusion_sort_label}"
        )
        return gold_standard_template, FP_genes


@dataclass
class CTAT(BaseFusion):
    def parse_results(self) -> pd.DataFrame:
        assert self.output_result_file.startswith("gs") or os.path.exists(
            self.output_result_file
        )
        self.result = pd.read_csv(self.output_result_file, sep="\t")
        print(self.result.head())
        return self.result

    @property
    def reads_field(self):
        return "num_LR"

    def clean_results(self):
        self.result["FusionName"] = self.result["#FusionName"].apply(
            lambda x: "--".join(map(lambda x: x.upper(), x.split("--")))
        )
        self.result["Breakpoint"] = self.result.loc[
            :, ["LeftBreakpoint", "RightBreakpoint"]
        ].apply(
            lambda x: "--".join(map(lambda y: ":".join(y.split(":")[:2]), x)), axis=1
        )
        self.result["FusionLexSort"] = self.result["#FusionName"].apply(
            lambda x: "--".join(sorted(map(lambda x: x.upper(), x.split("--"))))
        )
        self.result["BreakpointLexSort"] = self.result.loc[
            :, ["LeftBreakpoint", "RightBreakpoint"]
        ].apply(
            lambda x: "--".join(sorted(map(lambda y: ":".join(y.split(":")[:2]), x))),
            axis=1,
        )
        self.result = self.result.loc[
            self.result.loc[:, self.reads_field] > self.lower_threshold, :
        ]


@dataclass
class GMAPFusion(CTAT):
    @property
    def reads_field(self):
        return "JunctionReadCount"

    def parse_results(self) -> pd.DataFrame:
        assert os.path.exists(self.output_result_file)
        self.result = pd.read_table(self.output_result_file)
        return self.result


@dataclass
class Pbfusion(BaseFusion):
    def parse_results(self) -> pd.DataFrame:
        assert self.output_result_file.startswith("gs") or os.path.exists(
            self.output_result_file
        )

        # Different pbfusion version has different INFO lines
        # v0.3.0
        # INFO line number
        self.commented_rows = 20
        self.result = pd.read_table(
            self.output_result_file, skiprows=self.commented_rows
        )
        # except pd.errors.ParserError:
        if "info" not in self.result.columns:  # v0.1
            self.commented_rows = 18
            self.result = pd.read_table(
                self.output_result_file, skiprows=self.commented_rows
            )
        if "info" not in self.result.columns:  # v0.2.2
            self.commented_rows = 14
            self.result = pd.read_table(
                self.output_result_file, skiprows=self.commented_rows
            )
        assert "info" in self.result
        return self.result

    @property
    def reads_field(self):
        return "num_LR"

    def clean_results(self):
        print(self.output_result_file)
        pbfusion_read_count = self.result.loc[:, "info"].map(
            lambda x: int(x.split(";")[0].split("=")[-1])
        )
        self.result.loc[:, "num_LR"] = pbfusion_read_count

        if self.commented_rows in [18, 20]:
            sorted_gene_pairs = list(
                self.result.loc[:, "info"]
                .str.extract("(?<=GN=)([^;]+);")
                .iloc[:, 0]
                .str.upper()
                .map(lambda x: "--".join(sorted(x.split(","))))
            )
            unsorted_gene_pairs = list(
                self.result.loc[:, "info"]
                .str.extract("(?<=GN=)([^;]+);")
                .iloc[:, 0]
                .str.upper()
                .map(lambda x: "--".join(x.split(",")))
            )
        else:
            sorted_gene_pairs = list(
                self.result.loc[:, "info"]
                .str.extract("(?<=GENE_NAMES=)([^;]+);")
                .iloc[:, 0]
                .str.upper()
                .map(lambda x: "--".join(sorted(x.split(","))))
            )
            unsorted_gene_pairs = list(
                self.result.loc[:, "info"]
                .str.extract("(?<=GENE_NAMES=)([^;]+);")
                .iloc[:, 0]
                .str.upper()
                .map(lambda x: "--".join(x.split(",")))
            )

        self.result["FusionName"] = unsorted_gene_pairs
        self.result["FusionLexSort"] = sorted_gene_pairs

        self.result["Breakpoint"] = self.result.loc[
            :, ["#chr1", "start1", "chr2", "start2"]
        ].apply(lambda x: "--".join([f"{x[0]}:{x[1]}", f"{x[2]}:{x[3]}"]), axis=1)
        self.result["BreakpointLexSort"] = self.result.loc[
            :, ["#chr1", "start1", "chr2", "start2"]
        ].apply(
            lambda x: "--".join(sorted([f"{x[0]}:{x[1]}", f"{x[2]}:{x[3]}"])), axis=1
        )

        self.result = self.result.loc[pbfusion_read_count > self.lower_threshold, :]


@dataclass
class JAFFAL(BaseFusion):
    def parse_results(self) -> pd.DataFrame:
        assert self.output_result_file.startswith("gs") or os.path.exists(
            self.output_result_file
        )
        self.result = pd.read_csv(self.output_result_file)
        if self.jaffal_sample_id is not None:
            self.result = self.result.loc[
                self.result.iloc[:, 0].isin([self.jaffal_sample_id]), :
            ].copy()
        return self.result

    @property
    def reads_field(self):
        return "spanning reads"

    def clean_results(self):
        self.result["FusionName"] = (
            self.result["fusion genes"]
            .str.upper()
            .apply(lambda x: "--".join(x.split(":")))
        )
        self.result["Breakpoint"] = self.result.loc[
            :, ["chrom1", "base1", "chrom2", "base2"]
        ].apply(lambda x: "--".join([f"{x[0]}:{x[1]}", f"{x[2]}:{x[3]}"]), axis=1)
        self.result["FusionLexSort"] = (
            self.result["fusion genes"]
            .str.upper()
            .apply(lambda x: "--".join(sorted(x.split(":"))))
        )
        self.result["BreakpointLexSort"] = self.result.loc[
            :, ["chrom1", "base1", "chrom2", "base2"]
        ].apply(
            lambda x: "--".join(sorted([f"{x[0]}:{x[1]}", f"{x[2]}:{x[3]}"])), axis=1
        )
        self.result = self.result.loc[
            (self.result.loc[:, "spanning reads"] > self.lower_threshold), :
        ]


@dataclass
class FusionSeeker(BaseFusion):
    @property
    def reads_field(self):
        return "NumSupp"

    def parse_results(self) -> pd.DataFrame:
        assert self.output_result_file.startswith("gs") or os.path.exists(
            self.output_result_file
        )
        self.result = pd.read_table(self.output_result_file)
        return self.result

    def clean_results(self):
        self.result["FusionName"] = self.result.loc[:, ["Gene1", "Gene2"]].apply(
            lambda x: "--".join(x).upper(), axis=1
        )
        self.result["FusionLexSort"] = self.result["FusionName"].apply(
            lambda x: "--".join(sorted(x.split("--"))).upper()
        )
        self.result["Breakpoint"] = self.result.loc[
            :, ["Chrom1", "Breakpoint1", "Chrom2", "Breakpoint2"]
        ].apply(lambda x: "--".join([f"{x[0]}:{x[1]}", f"{x[2]}:{x[3]}"]), axis=1)
        self.result["BreakpointLexSort"] = self.result.loc[
            :, ["Chrom1", "Breakpoint1", "Chrom2", "Breakpoint2"]
        ].apply(
            lambda x: "--".join(sorted([f"{x[0]}:{x[1]}", f"{x[2]}:{x[3]}"])), axis=1
        )
        self.result = self.result.loc[
            (self.result.loc[:, self.reads_field] > self.lower_threshold), :
        ]

    # def parse_results(self):
    #     pass
    # @property
    # def reads_field(self):
    #     pass


@dataclass
class LongGF(BaseFusion):
    def parse_results(self) -> pd.DataFrame:
        assert self.output_result_file.startswith("gs") or os.path.exists(
            self.output_result_file
        )
        if self.output_result_file.startswith("gs"):
            if not os.path.exists(os.path.basename(self.output_result_file)):
                os.system(f"gsutil cp -r {self.output_result_file} .")
            self.output_result_file = os.path.basename(self.output_result_file)

        result_table = []
        with open(self.output_result_file) as fin:
            for line in fin:
                if "SumGF" in line:
                    elements = line.strip().split()
                    result_table.append(elements)
        self.result = pd.DataFrame(result_table)
        self.result.iloc[:, 2] = self.result.iloc[:, 2].astype(int)
        return self.result

    @property
    def reads_field(self):
        return "num_LR"

    def clean_results(self):
        self.result = self.result.drop([0], axis=1)
        self.result = self.result.rename(
            {1: "FusionName", 2: "num_LR", 3: "LeftBreakpoint", 4: "RightBreakpoint"},
            axis=1,
        )
        self.result["FusionName"] = (
            self.result["FusionName"]
            .str.upper()
            .apply(lambda x: "--".join(x.split(":")))
        )
        self.result["FusionLexSort"] = (
            self.result.iloc[:, 0]
            .str.upper()
            .map(lambda x: "--".join(sorted(x.split("--"))))
            .values
        )
        self.result.loc[:, "Breakpoint"] = (
            self.result.loc[:, ["LeftBreakpoint", "RightBreakpoint"]]
            .apply(lambda x: "--".join(x), axis=1)
            .values
        )
        self.result.loc[:, "BreakpointLexSort"] = (
            self.result.loc[:, ["LeftBreakpoint", "RightBreakpoint"]]
            .apply(lambda x: "--".join(sorted(x)), axis=1)
            .values
        )
        self.result = self.result.loc[
            (self.result.loc[:, self.reads_field] > self.lower_threshold), :
        ]


@dataclass
class FLAIRFUSION(BaseFusion):

    def parse_results(self) -> pd.DataFrame:
        assert self.output_result_file.startswith("gs") or os.path.exists(
            self.output_result_file
        )
        self.result = pd.read_table(self.output_result_file)
        return self.result

    @property
    def reads_field(self):
        return "spanning reads"

    def clean_results(self):
        self.result = self.result.rename(
            {"#name": "FusionName", "3' breakpoint": "LeftBreakpoint", "5' breakpoint": "RightBreakpoint"},
            axis=1,
        )
        self.result["FusionName"] = (
            self.result["FusionName"]
            .str.upper()
            .apply(lambda x: "--".join(x.split(":")))
        )
        self.result["FusionLexSort"] = (
            self.result.iloc[:, 0]
            .str.upper()
            .map(lambda x: "--".join(sorted(x.split("--"))))
            .values
        )
        #3'-PTPRT-chr20-42791192-0.604
        self.result.loc[:, "Breakpoint"] = (
            self.result.loc[:, ["LeftBreakpoint", "RightBreakpoint"]]
            .apply(lambda x: "--".join([":".join(x[0].split("-")[-3:-1]), ":".join(x[1].split("-")[-3:-1])]), axis=1)
            .values
        )
        print(self.result.head())
        self.result.loc[:, "BreakpointLexSort"] = (
            self.result.loc[:, "Breakpoint"]
            .apply(lambda x: "--".join(sorted(x.split("--"))))
            .values
        )
        print(self.result.head())
        self.result = self.result.loc[
            (self.result.loc[:, self.reads_field] > self.lower_threshold), :
        ]


def breakpoints_comparison(x, y):
    left1, right1 = x.split("--")
    left2, right2 = y.split("--")

    left1 = left1.split(":")
    left2 = left2.split(":")
    right1 = right1.split(":")
    right2 = right2.split(":")
    assert left1[0] == left2[0]
    assert right1[0] == right2[0]
    return np.abs(int(left1[1]) - int(left2[1])) + np.abs(
        int(right1[1]) - int(right2[1])
    )


def match_paralogs():
    """wget -c https://raw.githubusercontent.com/fusiontranscripts/FusionBenchmarking/b08aae832bc36dea4a874dad5bd230baa546fa18/resources/paralog_clusters.2020.I5.dat"""
    import urllib.request

    paralog_file = urllib.request.urlopen(
        "https://raw.githubusercontent.com/fusiontranscripts/FusionBenchmarking/b08aae832bc36dea4a874dad5bd230baa546fa18/resources/paralog_clusters.2020.I5.dat"
    )
    paralog_pairs_dict = defaultdict(list)
    paralogs = []
    n = 0
    for line in paralog_file:
        paralog_pair = line.strip().decode("utf-8").split()
        # NOTE: take uppper case to be consistent with gold standard and prediction
        paralog_pair = list(map(lambda x: x.upper(), paralog_pair))
        paralogs.append(paralog_pair)
        for paralog in paralog_pair:
            paralog_pairs_dict[paralog].append(n)
        n += 1

    for paralog in paralog_pairs_dict:
        assert len(paralog_pairs_dict[paralog]) != 0
    return paralog_pairs_dict, paralogs


def parse_breakpt(desc):
    desc = desc.split()
    one_gene_pair = re.search(r"(\S+)\|.*--(\S+)\|.*", desc[0])
    if one_gene_pair:
        one_gene_pair = list(one_gene_pair.groups())
    else:
        raise Exception("No gene pair matched")
    one_end = desc[2].split(":")
    another_end = desc[3].split(":")
    one_end_coord = re.search(r"(?<=-)(\S+)\[", one_end[1].split(",")[-1])
    if one_end_coord:
        one_end_coord = one_end_coord.groups()
    else:
        raise Exception
    return (
        one_end[0],
        int(one_end_coord[0]),
        another_end[0],
        int(another_end[1].split(",")[0].split("-")[0]),
        "--".join(one_gene_pair),
    )


def load_gold_standard_pbsim3(template_sequences):
    """
    Args
    --------
    :param template_sequences: `str`, fastq.gz paths
    """
    all_benchmark = dict()

    for folder in glob.glob(template_sequences):
        if ".csv" in folder:
            continue
        if ".txt" in folder:
            continue
        fastq = os.path.join(folder, "mix_ccs.fastq.gz")
        mixmap = os.path.join(folder, "mix_name_map.txt")
        print(folder, fastq, mixmap)

        read_names = []
        with gzip.open(fastq, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq-sanger"):
                desc = record.description
                read_names.append(desc.replace("/ccs", ""))
        print(read_names[:5])
        gold = pd.read_table(mixmap, header=None, index_col=0, sep="\t")
        gold = gold.reset_index()
        gold.loc[:, "read_names"] = gold.iloc[:, 0].str.replace(
            r"\/\d$", "", regex=True
        )
        gold = gold.drop_duplicates("read_names")
        gold.index = gold.read_names
        gold_counts = gold.loc[read_names, :].value_counts([1])
        gold_breakpoints_hg19 = (
            pd.DataFrame(gold_counts)
            .reset_index()
            .apply(lambda x: parse_breakpt(x[1]), axis=1)
            .tolist()
        )
        gold_breakpoints_hg19 = pd.DataFrame(gold_breakpoints_hg19)
        gold_breakpoints_hg38 = pd.concat(
            [
                gold_breakpoints_hg19.apply(
                    lambda x: ":".join(list(map(str, converter[x[0]][x[1]][0]))[:2]),
                    axis=1,
                ),
                gold_breakpoints_hg19.apply(
                    lambda x: ":".join(list(map(str, converter[x[2]][x[3]][0]))[:2]),
                    axis=1,
                ),
                gold_breakpoints_hg19.iloc[:, -1],
            ],
            axis=1,
        )
        gold_breakpoints_hg38.loc[:, "Breakpoint"] = gold_breakpoints_hg38.apply(
            lambda x: "--".join(x[:2]), axis=1
        )
        gold_breakpoints_hg38.loc[:, "BreakpointLexSort"] = gold_breakpoints_hg38.apply(
            lambda x: "--".join(sorted(x[:2])), axis=1
        )
        gold_breakpoints_hg38.loc[:, "FusionLexSort"] = gold_breakpoints_hg38.apply(
            lambda x: "--".join(sorted(x[4].upper().split("--"))), axis=1
        )
        gold_breakpoints_hg38 = gold_breakpoints_hg38.set_index(4)
        gold_breakpoints_hg38.insert(2, "FusionName", gold_breakpoints_hg38.index)
        gold_breakpoints_hg38.columns = [
            "LeftBreakPoint",
            "RightBreakPoint",
            "FusionName",
            "Breakpoint",
            "BreakpointLexSort",
            "FusionLexSort",
        ]
        print(gold_breakpoints_hg38.head())
        gold_breakpoints_hg38.loc[:, "num_LR"] = gold_counts.values
        # remove alternative contigs
        gold_breakpoints_hg38 = gold_breakpoints_hg38.loc[
            ~gold_breakpoints_hg38.BreakpointLexSort.str.contains("alt"), :
        ]
        gold_breakpoints_hg38 = gold_breakpoints_hg38.drop_duplicates(
            ["BreakpointLexSort", "FusionLexSort"]
        )

        gold_breakpoints_hg38 = gold_breakpoints_hg38.set_axis(
            gold_breakpoints_hg38.loc[:, "FusionLexSort"].str.upper(), axis="index"
        )
        gold_breakpoints_hg38.loc[:, "FusionName"] = gold_breakpoints_hg38.loc[
            :, "FusionName"
        ].str.upper()
        print(gold_breakpoints_hg38.head())
        all_benchmark[folder] = gold_breakpoints_hg38.copy()

    return all_benchmark


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


def load_reads(fastqs, gold_standard_breakpoints_hg38):
    """For each gene fusion pairs,
    return the gold standard read count by parsing the fastq.gz

    """
    all_benchmark = dict()
    all_benchmark_reads = dict()

    all_benchmark_with_breakpt = dict()
    for fastq in glob.glob(fastqs):
        sim_data_name = os.path.basename(fastq).replace(".fastq.gz", "")
        all_benchmark_with_breakpt[sim_data_name] = gold_standard_breakpoints_hg38.copy(
            deep=True
        )
        all_benchmark_with_breakpt[sim_data_name].loc[:, "num_LR"] = 0

        all_benchmark[sim_data_name] = defaultdict(int)
        all_benchmark_reads[sim_data_name] = defaultdict(list)
        with gzip.open(fastq, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq-sanger"):
                desc = record.description.split()
                one_gene_pair = re.search(r"(\S+)\|.*--(\S+)\|.*", desc[1])
                if one_gene_pair:
                    all_benchmark[sim_data_name][
                        "--".join(
                            sorted(map(lambda x: x.upper(), one_gene_pair.groups()))
                        )
                    ] += 1
                    all_benchmark_reads[sim_data_name][
                        "--".join(
                            sorted(map(lambda x: x.upper(), one_gene_pair.groups()))
                        )
                    ].append(desc[0])
                    all_benchmark_with_breakpt[sim_data_name].loc[
                        "--".join(
                            sorted(map(lambda x: x.upper(), one_gene_pair.groups()))
                        ),
                        "num_LR",
                    ] += 1
                else:
                    raise Exception

    # write tables of truth set fusion read counts
    if not os.path.exists("data"):
        os.makedirs("data")
    
    for sim_data_name in all_benchmark.keys():
        df_counts_dict = all_benchmark[sim_data_name]
        df = pd.DataFrame(
            list(df_counts_dict.items()), columns=["fusion_name", "num_reads"]
        )
        df["sim_data_name"] = sim_data_name
        df.to_csv(
            "data/" + sim_data_name + ".truthset_counts.tsv", sep="\t", index=False
        )

    # write tables of truth set fusion read counts
    for sim_data_name in all_benchmark_with_breakpt.keys():
        df = all_benchmark_with_breakpt[sim_data_name]
        df.loc[:, "read_names"] = "NA"
        reads = all_benchmark_reads[sim_data_name]
        for fusion_name in reads.keys():
            df.loc[fusion_name, "read_names"] = ";".join(reads[fusion_name])
        df["sim_data_name"] = sim_data_name
        df.to_csv(
            "data/" + sim_data_name + ".truthset_counts_with_breakpts.tsv",
            sep="\t",
            index=False,
        )
    return all_benchmark, all_benchmark_reads, all_benchmark_with_breakpt
