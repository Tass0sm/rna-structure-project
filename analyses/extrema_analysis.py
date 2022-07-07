import matplotlib.pyplot as plt
import pandas as pd
import plotting
import filtering
import averaging
import utils

###############################################################################
#                                     data                                    #
###############################################################################

# ctrl_dwell_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-dwell-times.csv")
# trizol_dwell_df = pd.read_csv("./data/trizol-old/trizol-dwell-times.csv")
# dna_dwell_df = pd.read_csv("./data/dna/dna_dwell_times.csv")

ctrl_phred_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-phred.csv")
# trizol_phred_df = pd.read_csv("./data/trizol-old/trizol_phred_df.csv")
# dna_phred_df = pd.read_csv("./data/dna/dna_phred_df.csv")
# annotations_df = pd.read_csv("./data/secondary_structure/annotations.csv")
labeling_df = pd.read_csv("./data/secondary_structure/hairpin_labeling.csv")

###############################################################################
#                               extrema analysis                              #
###############################################################################

def at_peak(diff_series, i):
    return diff_series.iloc[i - 1] > 0 and diff_series.iloc[i] < 0

def close_to_peak(diff_series, i):
    if at_valley(diff_series, i):
        return False

    window = range(i - 1, i + 2)
    return any(map(lambda j: at_peak(diff_series, j), window))

def at_valley(diff_series, i):
    return diff_series.iloc[i - 1] < 0 and diff_series.iloc[i] > 0

def close_to_valley(diff_series, i):
    if at_peak(diff_series, i):
        return False

    window = range(i - 1, i + 2)
    return any(map(lambda j: at_valley(diff_series, j), window))

def extrema_analysis():
    phred = ctrl_phred_df.mean().iloc[1:]
    phred_diff = phred.diff()

    label_curve = labeling_df.iloc[0, 2:]
    jerk_points = utils.series_jerk_points(label_curve)

    num_at_peak = sum(map(lambda i: close_to_peak(phred_diff, i), jerk_points))
    num_at_valley = sum(map(lambda i: close_to_valley(phred_diff, i), jerk_points))
    num_near_peak_and_valley = sum(map(lambda i: close_to_peak(phred_diff, i) and close_to_valley(phred_diff, i), jerk_points))

    # test: given a jerk point i, there a peak or valley in the phred curve within 1 base of i
    population_size = len(jerk_points)
    positive = num_at_peak + num_at_valley - num_near_peak_and_valley
    negative = population_size - positive

    print(f"Total population: {population_size}")
    print(f"Number of points which are near a peak XOR valley (exclusively): {positive}")
    print(f"Number of points which are not near a peak OR valley: {negative}")

    # other info
    total_of_peaks_and_valleys = sum(map(lambda i: at_peak(phred_diff, i) or at_valley(phred_diff, i), range(len(phred))))
    print(f"Total peaks and valleys in the phred curve: {total_of_peaks_and_valleys}")
    print(f"Percentage of peaks and valleys in the phred curve: {total_of_peaks_and_valleys / len(phred)}")

extrema_analysis()
