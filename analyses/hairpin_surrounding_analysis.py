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
#                          loop surrounding analysis                          #
###############################################################################

def loop_surrounding_analysis():
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

    phred = ctrl_phred_df.mean().iloc[1:]

    windows = averaging.get_windows_from_loops(phred, site_loops, 40, 60)
    curve = windows.mean() / phred.mean()
    curve.name = "ctrl9kb"

    plotting.plot_curves([curve], title="Window Analysis", output_file="new-func-test.png")

loop_surrounding_analysis()
