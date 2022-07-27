import matplotlib.pyplot as plt
import pandas as pd
import plotting
import filtering
import averaging
import parsing
import utils

###############################################################################
#                               extrema analysis                              #
###############################################################################

# utils

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

# peak -> bend analysis

def peak_to_bend_analysis():
    phred = ctrl_phred_df.mean().iloc[1:]
    phred_diff = phred.diff()

    peak_or_valley_points = list(filter(lambda i: at_peak(phred_diff, i) or at_valley(phred_diff, i), range(len(phred))))

    label_curve = labeling_df.iloc[0, 2:]
    jerk_points = utils.series_jerk_points(label_curve)

    num_bend_points = sum(map(lambda i: i in jerk_points, peak_or_valley_points))

    # test: given a jerk point i, there a peak or valley in the phred curve within 1 base of i
    population_size = len(peak_or_valley_points)
    positive = num_bend_points
    negative = population_size - positive

    print(f"Total peak or valley points: {population_size}")
    print(f"Number of points which are a bend: {positive}")
    print(f"Number of points which are not a bend: {negative}")

# bend -> peak analysis

def bend_to_peak_analysis(labeling_df, phred_df):
    phred = phred_df.mean().iloc[1:]
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

###############################################################################
#                            improving predictions                            #
###############################################################################

def improving_first_prediction():
    i = 0
    j = 60

    real_vienna = labeling.get_vienna("../my-structure-method-project/real_vienna.txt")[i:j+1]
    pred_vienna = labeling.get_vienna("../my-structure-method-project/pred_vienna.txt")[i:j+1]

    ctrl_phred_curve = plotting.get_curve(ctrl_phred_df, i, j, name="phred", method="average")
    pred_struct_curve = plotting.get_curve(pred_labeling_df, i, j, name="pred_struct", method="average")

    new_vienna = improving.improve_sequence(pred_vienna, pred_struct_curve, ctrl_phred_curve)
    new_pred_struct_curve = pd.Series(labeling.label_hairpins(new_vienna))
    new_pred_struct_curve.name = "improved"

    real_struct_curve = plotting.get_curve(real_labeling_df, i, j, name="real_struct", method="average")

    # pred_struct_curve = plotting.get_curve(pred_labeling_df, i, j, name="pred_struct", method="average")
    # ctrl_phred_curve = plotting.get_curve(ctrl_phred_df, i, j, name="phred", method="average")

    plotting.plot_curves([real_struct_curve,
                          pred_struct_curve,
                          new_pred_struct_curve,
                          ctrl_phred_curve],
                         title=f"Structure Prediction With Phred Enhancement",
                         output_file=f"./phred-prediction-enhancement.png")



###############################################################################
#                          loop surroundings analysis                         #
###############################################################################

# Manual

def selected_loop_surrounding_analysis(phred_df):
    site_loops = [(28, 34),
                  (73, 86),
                  (150, 157),
                  (260, 265),
                  (824, 827),
                  (1084, 1092),
                  (1228, 1233),
                  (1655, 1661),
                  (3860, 3868),
                  (3932, 3944),
                  (4462, 4468),
                  (4965, 4974),
                  (5506, 5511),
                  (7361, 7364),
                  (7449, 7452),
                  (7478, 7484),
                  (7865, 7871),
                  (7940, 7945),
                  (8000, 8005),
                  (8312, 8316),
                  (8585, 8592),
                  (8739, 8794),
                  (8922, 8928),
                  (9103, 9109)]

    phred = phred_df.mean().iloc[1:]

    windows = averaging.get_windows_from_loops(phred, site_loops, 40, 60)
    curve = windows.mean() / phred.mean()
    curve.name = "ctrl9kb"

    plotting.plot_curves([curve], title="Window Analysis", output_file="new-func-test.png")

# First Automatic Try

def automatic_loop_surrounding_analysis(labeling_df, phred_df):
    label_curve = labeling_df.iloc[0, 2:]
    phred = phred_df.mean().iloc[1:]

    patrick_loops = [(28, 34),
                     (73, 86),
                     (150, 157),
                     (260, 265),
                     (824, 827),
                     (1084, 1092),
                     (1228, 1233),
                     (1655, 1661),
                     (3860, 3868),
                     (3932, 3944),
                     (4462, 4468),
                     (4965, 4974),
                     (5506, 5511),
                     (7361, 7364),
                     (7449, 7452),
                     (7478, 7484),
                     (7865, 7871),
                     (7940, 7945),
                     (8000, 8005),
                     (8312, 8316),
                     (8585, 8592),
                     (8739, 8794),
                     (8922, 8928),
                     (9103, 9109)]

    patrick_windows = averaging.get_windows_from_loops(phred, patrick_loops, 30, 30)
    patrick_curve = patrick_windows.mean() / phred.mean()
    patrick_curve.name = "manual"

    site_loops = averaging.get_top_loops(label_curve)

    automatic_windows = averaging.get_windows_from_loops(phred, site_loops, 30, 30)
    automatic_curve = automatic_windows.mean() / phred.mean()
    automatic_curve.name = "automatic"

    plotting.plot_curves([automatic_curve,
                          patrick_curve], title="Window Analysis", output_file="comparison_figure.png")

# Second Automatic Try

# utils

# zero_patch_size = 4

# def parse_backward_until_zero_patch(s, i):
#     # unsafe
#     while not (label_curve.diff().iloc[i - zero_patch_size:i] == 0.0).all():
#         i -= 1

#         if i <= zero_patch_size:
#             break

#     return i

# def parse_until_zero_patch(s, i):
#     # unsafe
#     while not (label_curve.diff().iloc[i:i + zero_patch_size] == 0.0).all():
#         i += 1

#         if i >= len(s) - zero_patch_size:
#             break

#     return i

def get_stem_length(loop, label_curve):
    left_end = parsing.parse_backward_until_zero(label_curve, loop[0])
    right_end = parsing.parse_until_zero(label_curve, loop[1])

    left_length = loop[0] - left_end
    right_length = right_end - loop[1]

    return left_length, right_length

def stem_filter(n, label_curve):
    def f(loop):
        l = get_stem_length(loop, label_curve)

        return l[0] > n and l[1] > n

    return f

def second_automatic_loop_surrounding_analysis(labeling_df, phred_df):
    label_curve = labeling_df.iloc[0, 2:]
    phred = phred_df.mean().iloc[1:]

    patrick_loops = [(28, 34),
                     (73, 86),
                     (150, 157),
                     (260, 265),
                     (824, 827),
                     (1084, 1092),
                     (1228, 1233),
                     (1655, 1661),
                     (3860, 3868),
                     (3932, 3944),
                     (4462, 4468),
                     (4965, 4974),
                     (5506, 5511),
                     (7361, 7364),
                     (7449, 7452),
                     (7478, 7484),
                     (7865, 7871),
                     (7940, 7945),
                     (8000, 8005),
                     (8312, 8316),
                     (8585, 8592),
                     (8739, 8794),
                     (8922, 8928),
                     (9103, 9109)]

    patrick_windows = averaging.get_windows_from_loops(phred, patrick_loops, 30, 30)
    patrick_curve = patrick_windows.mean() / phred.mean()
    patrick_curve.name = "manual"

    site_loops = averaging.get_top_loops(label_curve)

    automatic_windows = averaging.get_windows_from_loops(phred, site_loops, 30, 30)
    automatic_curve = automatic_windows.mean() / phred.mean()
    automatic_curve.name = "automatic"

    filtered_loops = list(filter(stem_filter(10, label_curve), site_loops))

    filtered_windows = averaging.get_windows_from_loops(phred, filtered_loops, 30, 30)
    filtered_curve = filtered_windows.mean() / phred.mean()
    filtered_curve.name = "filtered"

    plotting.plot_curves([automatic_curve,
                          filtered_curve,
                          patrick_curve], title="Window Analysis", output_file="comparison_figure.png")

# Third Automati Try

def third_automatic_loop_surrounding_analysis(labeling_df, phred_df):
    label_curve = labeling_df.iloc[0, 2:]
    phred = phred_df.mean().iloc[1:]

    patrick_loops = [(28, 34),
                     (73, 86),
                     (150, 157),
                     (260, 265),
                     (824, 827),
                     (1084, 1092),
                     (1228, 1233),
                     (1655, 1661),
                     (3860, 3868),
                     (3932, 3944),
                     (4462, 4468),
                     (4965, 4974),
                     (5506, 5511),
                     (7361, 7364),
                     (7449, 7452),
                     (7478, 7484),
                     (7865, 7871),
                     (7940, 7945),
                     (8000, 8005),
                     (8312, 8316),
                     (8585, 8592),
                     (8739, 8794),
                     (8922, 8928),
                     (9103, 9109)]

    patrick_windows = averaging.get_windows_from_loops(phred, patrick_loops, 30, 30)
    patrick_curve = patrick_windows.mean() / phred.mean()
    patrick_curve.name = "manual"

    site_loops = averaging.get_top_loops(label_curve)

    automatic_windows = averaging.get_windows_from_loops(phred, site_loops, 30, 30)
    automatic_curve = automatic_windows.mean() / phred.mean()
    automatic_curve.name = "automatic"

    filtered_loops = list(filter(stem_filter(10, label_curve), site_loops))

    filtered_windows = averaging.get_windows_from_loops(phred, filtered_loops, 30, 30)
    filtered_curve = filtered_windows.mean() / phred.mean()
    filtered_curve.name = "filtered"

    plotting.plot_curves([automatic_curve,
                          filtered_curve,
                          patrick_curve], title="Window Analysis", output_file="comparison_figure.png")


###########################################################################
#                                Debugging                                #
###########################################################################

def loop_result_comparison(labeling_df, patrick_df):
    label_curve = labeling_df.iloc[0, 2:]

    patrick_loops = [(28, 34),
                     (73, 86),
                     (150, 157),
                     (260, 265),
                     (824, 827),
                     (1084, 1092),
                     (1228, 1233),
                     (1655, 1661),
                     (3860, 3868),
                     (3932, 3944),
                     (4462, 4468),
                     (4965, 4974),
                     (5506, 5511),
                     (7361, 7364),
                     (7449, 7452),
                     (7478, 7484),
                     (7865, 7871),
                     (7940, 7945),
                     (8000, 8005),
                     (8312, 8316),
                     (8585, 8592),
                     (8739, 8794),
                     (8922, 8928),
                     (9103, 9109)]
    site_loops = averaging.get_top_loops(label_curve)
    filtered_loops = list(filter(stem_filter(10, label_curve), site_loops))

    problem_loops = set(filtered_loops).difference(set(patrick_loops))

    for loop in problem_loops:
        print(loop)
        i = max(loop[0] - 100, 0)
        j = min(loop[1] + 100, 9173)
        label_curve = plotting.get_curve(labeling_df, i, j, name="struct", method="average")
        annotation_curve = plotting.get_curve(patrick_df, i, j, name="patrick", method="average")

        plotting.plot_curves([label_curve,
                              annotation_curve],
                             title=f"Structure around Problem Loop {i}-{j}",
                             output_file=f"./{i}-{j}.png")

def unfound_loop_analysis(labeling_df, patrick_df):
    label_curve = labeling_df.iloc[0, 2:]

    patrick_loops = [(28, 34),
                     (73, 86),
                     (150, 157),
                     (260, 265),
                     (824, 827),
                     (1084, 1092),
                     (1228, 1233),
                     (1655, 1661),
                     (3860, 3868),
                     (3932, 3944),
                     (4462, 4468),
                     (4965, 4974),
                     (5506, 5511),
                     (7361, 7364),
                     (7449, 7452),
                     (7478, 7484),
                     (7865, 7871),
                     (7940, 7945),
                     (8000, 8005),
                     (8312, 8316),
                     (8585, 8592),
                     (8739, 8794),
                     (8922, 8928),
                     (9103, 9109)]
    site_loops = averaging.get_top_loops(label_curve)
    filtered_loops = list(filter(stem_filter(10, label_curve), site_loops))

    problem_loops = set(patrick_loops).difference(set(filtered_loops))

    for loop in problem_loops:
        print(loop)
        i = max(loop[0] - 100, 0)
        j = min(loop[1] + 100, 9173)
        label_curve = plotting.get_curve(labeling_df, i, j, name="struct", method="average")
        annotation_curve = plotting.get_curve(patrick_df, i, j, name="patrick", method="average")

        plotting.plot_curves([label_curve,
                              annotation_curve],
                             title=f"Structure around Problem Loop {i}-{j}",
                             output_file=f"./{i}-{j}.png")

def loop_surrounding_debugging(labeling_df):
    label_curve = labeling_df.iloc[0, 2:]
    site_loops = averaging.get_top_loops(label_curve)

    for loop in site_loops:
        i = loop[0] - 20
        j = loop[1] + 20
        curve = plotting.get_curve(labeling_df, i, j, name="struct", method="average")

        plotting.plot_curves([curve],
                             title=f"Phred Score Curve Around Loop {i}-{j}",
                             output_file=f"./{i}-{j}.png")
