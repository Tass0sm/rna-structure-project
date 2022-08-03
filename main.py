import matplotlib.pyplot as plt
import pandas as pd
import loading
import plotting
import parsing
import labeling
import filtering
import analyses
import averaging
import utils

import improving

###############################################################################
#                                     data                                    #
###############################################################################

# ctrl_dwell_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-dwell-times.csv")
# trizol_dwell_df = pd.read_csv("./data/trizol-old/trizol-dwell-times.csv")
# dna_dwell_df = pd.read_csv("./data/dna/dna_dwell_times.csv")

seq_df = pd.read_csv("./data/sequence.csv")
# ctrl_phred_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-phred.csv")
# trizol_phred_df = pd.read_csv("./data/trizol-old/trizol_phred_df.csv")
# dna_phred_df = pd.read_csv("./data/dna/dna_phred_df.csv")
# patrick_df = pd.read_csv("./data/secondary_structure/annotations.csv")
# real_labeling_df = pd.read_csv("./data/secondary_structure/hairpin_labeling.csv")
# pred_labeling_df = pd.read_csv("./data/secondary_structure/pred_labeling.csv")

###############################################################################
#                                     old                                     #
###############################################################################

seq_curve = seq_df.iloc[:, 0]
cols = [seq_curve]

for n in ["./data/false_peaks/NEW/newF1F2-avg-phred.csv",
          "./data/false_peaks/NEW/avg-dwell-times.csv",
          "./data/false_peaks/OLD/OLD_IVT-avg-phred.csv",
          "./data/false_peaks/OLD/avg-dwell-times.csv",
          "./data/ctrl9kb/ctrl9kb-avg-phred.csv",
          "./data/ctrl9kb/avg-dwell-times.csv"]:
    df = pd.read_csv(n)
    col_s = df.iloc[:, 1]
    col_s.name = n
    cols.append(col_s)

false_peaks_result_df = pd.concat(cols, axis=1)
false_peaks_result_df.to_excel("./false_peaks_result.xlsx")

# # First

# mutation8079_df = pd.read_csv("./data/structure_interaction/8079mutation_phred.csv")
# WTcellular_8079MOD_df = pd.read_csv("./data/structure_interaction/WTcellular_8079MOD_phred.csv")
# WTcellular_8079UNM_df = pd.read_csv("./data/structure_interaction/WTcellular_8079UNM_phred.csv")

# p = 8042
# i = p - 20
# j = p + 20

# mutation8079_curve = plotting.get_curve(mutation8079_df, i, j, name="mutation8079", method="average")
# WTcellular_8079MOD_curve = plotting.get_curve(WTcellular_8079MOD_df, i, j, name="WTcellular_8079MOD", method="average")
# WTcellular_8079UNM_curve = plotting.get_curve(WTcellular_8079UNM_df, i, j, name="WTcellular_8079UNM", method="average")

# plotting.plot_curves([mutation8079_curve,
#                       WTcellular_8079MOD_curve,
#                       WTcellular_8079UNM_curve],
#                      title=f"Phred score around {p}",
#                      output_file=f"./phred-change-around-{p}.png")

# # Second

# mutation8989_df = pd.read_csv("./data/structure_interaction/8989mutation_phred.csv")
# WTcellular_8989MOD_df = pd.read_csv("./data/structure_interaction/WTcellular_8989MOD_phred.csv")
# WTcellular_8989UNM_df = pd.read_csv("./data/structure_interaction/WTcellular_8989UNM_phred.csv")

# p = 8797
# i = p - 20
# j = p + 20

# mutation8989_curve = plotting.get_curve(mutation8989_df, i, j, name="mutation8989", method="average")
# WTcellular_8989MOD_curve = plotting.get_curve(WTcellular_8989MOD_df, i, j, name="WTcellular_8989MOD", method="average")
# WTcellular_8989UNM_curve = plotting.get_curve(WTcellular_8989UNM_df, i, j, name="WTcellular_8989UNM", method="average")

# plotting.plot_curves([mutation8989_curve,
#                       WTcellular_8989MOD_curve,
#                       WTcellular_8989UNM_curve],
#                      title=f"Phred score around {p}",
#                      output_file=f"./phred-change-around-{p}.png")

# i = 0
# j = 9173

# mutation8079_df = pd.read_csv("./data/structure_interaction/8079mutation_phred.csv")
# WTcellular_8079MOD_df = pd.read_csv("./data/structure_interaction/WTcellular_8079MOD_phred.csv")
# WTcellular_8079UNM_df = pd.read_csv("./data/structure_interaction/WTcellular_8079UNM_phred.csv")
# mutation8989_df = pd.read_csv("./data/structure_interaction/8989mutation_phred.csv")
# WTcellular_8989MOD_df = pd.read_csv("./data/structure_interaction/WTcellular_8989MOD_phred.csv")
# WTcellular_8989UNM_df = pd.read_csv("./data/structure_interaction/WTcellular_8989UNM_phred.csv")

# mutation8079_curve = plotting.get_curve(mutation8079_df, i, j, name="mutation8079", method="average")
# WTcellular_8079MOD_curve = plotting.get_curve(WTcellular_8079MOD_df, i, j, name="WTcellular_8079MOD", method="average")
# WTcellular_8079UNM_curve = plotting.get_curve(WTcellular_8079UNM_df, i, j, name="WTcellular_8079UNM", method="average")
# mutation8989_curve = plotting.get_curve(mutation8989_df, i, j, name="mutation8989", method="average")
# WTcellular_8989MOD_curve = plotting.get_curve(WTcellular_8989MOD_df, i, j, name="WTcellular_8989MOD", method="average")
# WTcellular_8989UNM_curve = plotting.get_curve(WTcellular_8989UNM_df, i, j, name="WTcellular_8989UNM", method="average")


# curves = [seq_curve, WTcellular_8989UNM_curve, WTcellular_8989MOD_curve, mutation8989_curve, WTcellular_8079UNM_curve, WTcellular_8079MOD_curve, mutation8079_curve]
