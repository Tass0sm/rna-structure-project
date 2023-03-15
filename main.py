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

# Fast5s: /fs/project/PAS1405/General/HIV_RNA_modification_dataset/HIVDNA/072517_HIVDNA_reverse_filtered_single_fast5/
# fastq: /fs/project/PAS1405/General/HIV_RNA_modification_dataset/HIVDNA/072517_HIVDNA_reverse.fastq

# rna_seq_df = pd.read_csv("./data/rna_seq_df.csv")
# dna_seq_df = pd.read_csv("./data/dna_seq_df.csv")
# ctrl_phred_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-phred.csv")
# trizol_phred_df = pd.read_csv("./data/trizol-old/trizol_phred_df.csv")
# dna_phred_df = pd.read_csv("./data/dna_sense/dna_sense_phreds.csv")
dd_pcr_forward_phred_df = pd.read_csv("./data/dd_pcr_forward/dd_pcr_forward_phreds.csv")
dd_pcr_reverse_phred_df = pd.read_csv("./data/dd_pcr_reverse/dd_pcr_reverse_phreds.csv")
# reverse_dna_phred_df = pd.read_csv("./data/dna/dna-alignment-with-rna-aligned-phreds.csv")
# patrick_df = pd.read_csv("./data/secondary_structure/annotations.csv")
# real_labeling_df = pd.read_csv("./data/secondary_structure/hairpin_labeling.csv")
# pred_labeling_df = pd.read_csv("./data/secondary_structure/pred_labeling.csv")

###############################################################################
#                                     old                                     #
###############################################################################

def loop_surrounding_analysis(phred_df):
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
    return patrick_windows


# dd_pcr_reverse_phred_df

output_dir = "./data/dd_pcr_forward/"
base_name = "dd_pcr_forward"

patrick_windows = loop_surrounding_analysis(dd_pcr_forward_phred_df)
patrick_windows.to_excel(output_dir + base_name + "_phreds_around_patrick_loops.xlsx")

dna_seq_df = pd.read_csv("./data/dna_seq_df.csv")
dna_seq_df = dna_seq_df.set_index("Unnamed: 0")
mean_phred_df = pd.DataFrame(dd_pcr_forward_phred_df.mean().iloc[1:]).T
combined_df = pd.concat([dna_seq_df, mean_phred_df], axis=0).T
combined_df.to_excel(output_dir + base_name + "_average_phreds_with_sequence.xlsx")

# def loop_surrounding_analysis(phred_df):
#     phred = phred_df.mean().iloc[1:]



#     patrick_windows = averaging.get_windows_from_loops(phred, patrick_loops, 30, 30)
#     return patrick_windows




    # patrick_curve = patrick_windows.mean() / phred.mean()
    # patrick_curve.name = "manual"

    # plotting.plot_curves([automatic_curve,
    #                       patrick_curve], title="Window Analysis", output_file="comparison_figure.png")


# for n in ["IVT/newF1F2_GL",
#           "OLD_IVT/OLD_IVT"]:
#     #f = f"/fs/project/PAS1405/General/HIV_RNA_modification_dataset/{n}.fastq"
#     o = f"/users/PAS1405/tassosm/Desktop/rna-structure-project/data/false_peaks/{n}.sam"
#     #loading.align_fastq(f, o)
#     df = loading.read_aligned_sequences(o, limit=8192)
#     c = f"/users/PAS1405/tassosm/Desktop/rna-structure-project/data/false_peaks/{n}.csv"
#     df.to_csv(c)

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
