import pandas as pd

from sys import argv
import scipy.ndimage.filters as filters
import scipy.signal as sig

from sklearn.utils.fast_dict import IntFloatDict as fdict

def create_fastdict(keys=None, vals=None):
    if keys is None:
        d = fdict(np.array([], dtype=int), np.array([]))
        return d
    else:
        return fdict(keys, vals)

import math
sqrt = math.sqrt
import numpy as np

from time import time, sleep


def _apply_filter(arr, width):

    "should probably memoize this function :"

    val = (2 * width + 1)
    moving_avg_a = np.ones(val) * 1/(val)
    # start = time()
    # print("get_window")
    a = sig.get_window('hanning', val)

    # print("convolve1d")
    ga = filters.convolve1d(arr, a/a.sum())

    # print("argrelextrema")
    minima_a = sig.argrelextrema(ga, np.less)[0]

    return minima_a


def _find_end(a, init_search_location, n_bpoints):
    search_location = init_search_location
    filtered = _apply_filter(a, search_location)
    while len(filtered) >= n_bpoints:
        filtered = _apply_filter(a, search_location * 2)
        search_location = search_location * 2

    return search_location


def _width_raw(a, n_bpoints):

    "gives the same result as the variable found_width in the original code"

    search_val = n_bpoints

    filtered = _apply_filter(a, search_val)
    idx = np.searchsorted(filtered, search_val, side="left") - 1

    return len(filtered) - idx


def _trackback(array, search_val, start_search):

    step_coarse = 20
    delta_coarse = 200
    step_fine = 1

    for i in range(start_search + step_coarse, len(array), step_coarse):

        filtered = _apply_filter(array, i)

        if len(filtered) == search_val:
            start_search = i


    for i in range(start_search + step_fine, len(array), step_fine):

        filtered = _apply_filter(array, i)

        if len(filtered) == search_val:
            start_search = i

    return start_search


def _width(val, n_bpoints):

    search_value = 1000
    end = _find_end(val, search_value, n_bpoints)

    found_width = end - _width_raw(val, n_bpoints)

    found_width_trackback_raw = end - _trackback(val, search_value, found_width)

    return found_width_trackback_raw


def find_breakpoint_loci(df, n_snps=50):

    n_bpoints = int(math.ceil(len(df) / n_snps) - 1)

    val = df.val
    pos = df.pos.values

    width = _width(val, n_bpoints)

    filtered = _apply_filter(val, width)
    return pos[filtered]



from helpers import (covar_files_map, preread_files, find_start_locus, find_end_locus,
                     update_covar_and_loci)


def compute_metrics(covariance_files, partitions, breakpoints):

    snp_first, snp_last = partitions.head(1).Start.iloc[0], partitions.tail(1).End.iloc[0]

    curr_locus = -1

    covar_dict = covar_files_map(covariance_files, partitions)

    covar, loci, iter_start = preread_files(partitions, covar_dict, snp_first)

    curr_locus = -1
    end_locus = 0

    curr_breakpoint_index = 0
    curr_breakpoint = breakpoints[curr_breakpoint_index]

    curr_locus_index = 0
    block_height = 0
    block_width = 0
    block_width_sum = 0
    total_snps = 0
    nzero = 0
    nonzero = 0
    metric_sum = 0


    for i in range(iter_start, len(partitions)):

        j = 0

        pstart, pend = partitions.iloc[i]
        f = covar_dict[pstart, pend]
        covar, loci, covar_starts, covar_ends = update_covar_and_loci(covar, loci, end_locus, f)

        autocovar = covar[covar.i.values == covar.j.values].drop("j", 1).set_index("i").squeeze().to_dict()

        i_, j_, vals = [covar[v].values for v in "i j val".split()]

        covar_len = len(covar)

        start_locus_index, start_locus = find_start_locus(curr_locus, loci, snp_first)
        curr_locus_index, curr_locus = start_locus_index, start_locus
        end_locus_index, end_locus = find_end_locus(i, partitions, loci, snp_last)

        while curr_locus <= end_locus:
            if curr_locus > curr_breakpoint:

                block_height = 0 - total_snps
                nzero += block_height * block_width
                block_width_sum += block_width

                block_width = 0
                curr_breakpoint_index += 1

                if curr_breakpoint_index >= len(breakpoints):
                    break

                curr_breakpoint = breakpoints[curr_breakpoint_index]

            block_width += 1

            # anything with j is just more convoluted code to get the covars more quickly
            while j < covar_len and i_[j] != curr_locus:
                j += 1

            while j < covar_len and i_[j] == curr_locus:
                if j_[j] > curr_breakpoint:
                    corrcoeff = vals[j] / sqrt(autocovar[i_[j]] * autocovar[j_[j]])
                    metric_sum += corrcoeff ** 2
                    nonzero += 1

                j += 1

            if curr_locus_index + 1 < len(loci):
                curr_locus_index += 1
                curr_locus = loci[curr_locus_index]
                total_snps += 1
            else:
                break

    nzero += (total_snps * block_width_sum)
    # df = pd.DataFrame({"Names": ["metric_sum", "nzero", "nonzero"], "Values": [metric_sum, nzero, nonzero]})
    # m = {}
    # m["metric_sum"]

    return metric_sum, nzero, nonzero

if __name__ == "__main__":

    f = argv[1]

    df = pd.read_table(f, header=None, names="pos val".split())

    partition_file = argv[2]

    partitions = pd.read_table(partition_file, sep=" ", header=None)
    partitions.index = range(len(partitions))
    partitions.columns = ["Start", "End"]


    breakpoints = find_breakpoint_loci(df)

    covariance_files = argv[3:]

    metrics = compute_metrics(covariance_files, partitions, breakpoints)

    # df.to_csv("metrics.tsv", sep="\t", index=False)
    # pd.Series(breakpoints).to_frame().to_csv("breakpoints.tsv", sep="\t", index=False, header=False)
