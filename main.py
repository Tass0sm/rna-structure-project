import matplotlib.pyplot as plt
import pandas as pd
import plotting
import filtering
import averaging

labeling_df = pd.read_csv("./data/secondary_structure/hairpin_labeling.csv")
structure_curve = labeling_df.iloc[0, 2:]
structure_diff = structure_curve.diff()
# jerk_points = plotting.series_jerk(label_curve)
# jerk_points = jerk_points[jerk_points != 0]
# points = list(map(lambda x: int(x), jerk_points.index))

ctrl_phred_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-phred.csv")
phred = ctrl_phred_df.mean().iloc[1:]

# ctrl_dwell_df = pd.read_csv("./data/ctrl9kb/ctrl9kb-dwell-times.csv")
# trizol_phred_df = pd.read_csv("./data/trizol-old/trizol_phred_df.csv")
# trizol_dwell_df = pd.read_csv("./data/trizol-old/trizol-dwell-times.csv")
# dna_phred_df = pd.read_csv("./data/dna/dna_phred_df.csv")
# dna_dwell_df = pd.read_csv("./data/dna/dna_dwell_times.csv")

secondary_structure_df = pd.read_csv("./data/secondary_structure/secondary_structure.csv")
struct = secondary_structure_df.iloc[0, 2:]
sites = plotting.get_sites(struct)

def get_site_loop(i):
    _, loop = averaging.parse_loop(structure_diff, i)
    return loop

site_loops = list(map(lambda s: get_site_loop(s[0]), sites))

windows = averaging.get_windows_from_loops(phred, site_loops, 30, 30)
curve = windows.mean().values / phred.mean()

fig, ax = plt.subplots(figsize=(12,5))

# for p in points[0:10]:
#     # df = averaging.get_windows(sequence, points, -10, 10)
#     window = averaging.safe_get_window(sequence, p - 10, p + 10)

#     fig, ax = plt.subplots(figsize=(12,5))

ax.plot(list(range(-30, 31)), curve)

# ax.set_xticks(np.arange(beginning_pos, ending_pos, 25))
# ax.set_xticklabels(np.arange(beginning_pos, ending_pos, 25))
ax.set_xlabel('position (bp)')
# ax.set_xlim((beginning_pos, ending_pos))
# ax.set_ylim((0, 40))

ax.set_title("Dr. Kim Analysis")
fig.savefig(f"./dr-kim-analysis.png")
plt.close(fig)

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
