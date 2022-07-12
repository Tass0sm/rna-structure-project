import matplotlib.pyplot as plt
import pandas as pd
import plotting
import filtering
import analyses
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
#                                     old                                     #
###############################################################################

analyses.second_automatic_loop_surrounding_analysis(labeling_df, ctrl_phred_df)

# struct = secondary_structure_df.iloc[0, 2:]
# sites = plotting.get_sites(struct)

# def get_site_loop(i):
#     _, loop = averaging.parse_loop(structure_diff, i)
#     return loop

# site_loops = averaging.get_loops(structure_curve.diff())
# site_loops = list(map(lambda s: get_site_loop(s[0]), sites))

# site_loops = [(28, 34),
#               (73, 86),
#               (150, 157),
#               (254, 263),
#               (824, 827),
#               (1084, 1092),
#               (1201, 1204),
#               (1617, 1622),
#               (3860, 3904),
#               (4460, 4470),
#               (4965, 4974),
#               (5506, 5511),
#               (7449, 7452),
#               (7476, 7484),
#               (7865, 7880),
#               (8733, 8740)]



# sequence = pd.read_csv("./data/sequence.csv")

# avg_phred_per_kmer_df = pd.read_csv("/home/tassos/work/old/structure-data-project/cell/avg_phred_per_kmer.csv")
# avg_phred_per_kmer_df.set_index(avg_phred_per_kmer_df.columns[0], inplace=True)
# kmer_dict = dict(zip(avg_phred_per_kmer_df.index.values, avg_phred_per_kmer_df.iloc[:, 0].values))

# for i, j in sites:
#     ctrl_phred_curve = plotting.get_curve(ctrl_phred_df, i, j, name="ctrl9kb", method="average")
#     # ctrl_phred_stddev_curve = plotting.get_curve(ctrl_phred_df, i, j, name="ctrl9kb-stddev", method="stddev")

#     # normalized_ctrl_phred_curve = filtering.normalize_by_kmer(ctrl_phred_curve, sequence, kmer_dict)
#     # normalized_ctrl_phred_curve.name = "normalized"
#     structure_curve = plotting.get_curve(secondary_structure_df, i, j, name="struct", method="average")
#     label_curve = plotting.get_curve(labeling_df, i, j, name="struct", method="average")

#     plotting.custom_plot_curves([ctrl_phred_curve,
#                                  structure_curve,
#                                  label_curve],

#                                 [plotting.add_curve,
#                                  plotting.add_curve,
#                                  plotting.add_label_curve],

#                                 title="Phred Score Curve Aligned to Secondary Structure Numerical Representation",
#                                 output_file=f"plots/{i}-{j}.png")

# TAR-PBS site (positions 1 to 430)
# See Fig.3 https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004230
# plot_curves([("ctrl9kb", ctrl_phred_df)],
#             1, 430,
#             method="average",
#             title="TAR-PBS - ctrl9kb Average Phred per Base",
#             output_file="test-TAR-PBS.png")

# # TAR-PolA (positions 9000-9173)
# # See Fig.4 https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004230
# plot_curves([("ctrl9kb", ctrl_phred_df),
#              ("trizol", trizol_phred_df),
#              ("dna", dna_phred_df),
#              ("struct", secondary_structure_df)],
#             9000, 9173,
#             average=True,
#             title="TAR-PolA - ctrl9kb, trizol, dna, and secondary structure",
#             output_file="structure-TAR-PolA.png")

# # RRE (rev responsive element) (positions 7101-7700)
# # See Fig.4
# plot_curves([("ctrl9kb", ctrl_phred_df),
#              ("trizol", trizol_phred_df),
#              ("dna", dna_phred_df),
#              ("struct", secondary_structure_df)],
#             7101, 7700,
#             average=True,
#             title="RRE - ctrl9kb, trizol, dna, and secondary structure",
#             output_file="structure-RRE.png")

# # near PPT (positions 8402-8805)
# # See Fig.4
# plot_curves([("ctrl9kb", ctrl_phred_df),
#              ("trizol", trizol_phred_df),
#              ("dna", dna_phred_df),
#              ("struct", secondary_structure_df)],
#             8402, 8805,
#             average=True,
#             title="near-PPT - ctrl9kb, trizol, dna, and secondary structure",
#             output_file="structure-near-PPT.png")

# # near 8000 (positions 7700-8385)
# # See Fig.4
# plot_curves([("ctrl9kb", ctrl_phred_df),
#              ("trizol", trizol_phred_df),
#              ("dna", dna_phred_df),
#              ("struct", secondary_structure_df)],
#             7700, 8385,
#             average=True,
#             title="near-8000 - ctrl9kb, trizol, dna, and secondary structure",
#             output_file="structure-near-8000.png")

# # near 8988 (positions 8988-9020)
# # See Fig.4
# plot_curves([("ctrl9kb", ctrl_phred_df),
#              ("trizol", trizol_phred_df),
#              ("dna", dna_phred_df),
#              ("struct", secondary_structure_df)],
#             8802, 9020,
#             average=True,
#             title="near-8988 - ctrl9kb, trizol, dna, and secondary structure",
#             output_file="structure-near-8988.png")
