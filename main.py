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

# seq_df = pd.read_csv("./data/sequence.csv")
# ctrl_phred_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-phred.csv")
# trizol_phred_df = pd.read_csv("./data/trizol-old/trizol_phred_df.csv")
# dna_phred_df = pd.read_csv("./data/dna/dna_phred_df.csv")
# patrick_df = pd.read_csv("./data/secondary_structure/annotations.csv")
# real_labeling_df = pd.read_csv("./data/secondary_structure/hairpin_labeling.csv")
# pred_labeling_df = pd.read_csv("./data/secondary_structure/pred_labeling.csv")

###############################################################################
#                                     old                                     #
###############################################################################




for n in ["NEW/newF1F2_GL",
          "OLD/OLD_IVT"]:
    # f = f"./data/false_peaks/{n}.fastq"
    o = f"./data/false_peaks/{n}.sam"
    # loading.align_fastq(f, o)
    df = loading.read_aligned_sequences(o)
    c = f"./data/false_peaks/{n}.csv"
    df.to_csv(c)

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

# seq_curve = seq_df.iloc[:, 0]

# curves = [seq_curve, WTcellular_8989UNM_curve, WTcellular_8989MOD_curve, mutation8989_curve, WTcellular_8079UNM_curve, WTcellular_8079MOD_curve, mutation8079_curve]
# avg_phred_df = pd.concat(curves, axis=1)
# avg_phred_df.to_excel("./average_aligned_phred_scores.xlsx")
