import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from Bio import SeqIO

from itertools import product
from pysam import AlignmentFile
from tqdm import tqdm

###############################################################################
#                             Getting Phred Scores                            #
###############################################################################

def read_aligned_sequences(filename):
    alignment = AlignmentFile(filename)
    ref_len = alignment.lengths[0]
    read_len = alignment.count()
    alignment.close()

    col_index = ["read_id"] + list(range(ref_len))

    results_df = pd.DataFrame(index=range(read_len),
                              columns=col_index)

    alignment = AlignmentFile(filename)

    for i, record in tqdm(enumerate(alignment.fetch())):
        results_df.iloc[i, 0] = record.qname

        if (record.qual is not None):
            for p, q in zip(record.get_reference_positions(full_length=True),
                            record.qual):
                if (p is not None):
                    phred = ord(q) - 33
                    results_df.iloc[i, p + 1] = phred

    alignment.close()

    return results_df

def get_curve(df, begin, end, name="untitled", method="average"):
    subset_df = df.iloc[:, (begin + 2):(end + 3)]
    # print(f"Number of datapoints: {df.count().sum()}")

    if method == "average":
        curve_df = subset_df.mean()
    elif method == "minimum":
        curve_df = subset_df.min()
    elif method == "stddev":
        curve_df = subset_df.std()
    else:
        curve_df = subset_df.T

    curve_df.name = name
    beginning_pos = int(curve_df.index[0])
    ending_pos = int(curve_df.index[-1])
    curve_df.index = pd.RangeIndex(start=beginning_pos, stop=ending_pos + 1)

    return curve_df

def add_curve(ax, curve):
    curve.plot(ax = ax, legend = True, label = curve.name)

def plot_curves(curves, title="Per Base Sequence Quality", output_file="aligned_plot.png"):
    beginning_pos = int(curves[0].index[0])
    ending_pos = int(curves[0].index[-1])
    length = ending_pos - beginning_pos
    fig, ax = plt.subplots(figsize=(12,5))

    # Adding the lines.
    for curve in curves:
        add_curve(ax, curve)

    ax.set_xticks(np.arange(beginning_pos, ending_pos, 25))
    ax.set_xticklabels(np.arange(beginning_pos, ending_pos, 25))
    ax.set_xlabel('position (bp)')
    ax.set_xlim((beginning_pos, ending_pos))
    ax.set_ylim((-20, 40))

    ax.set_title(title)
    fig.savefig(output_file)
    plt.close(fig)

###############################################################################
#                           Normalized phred curves                           #
###############################################################################

def get_seq_dict(filename):
    global seq_dict

    seq_dict = {}
    alignment = AlignmentFile(filename)

    for record in tqdm(alignment.fetch()):
        seq_dict[record.qname] = record.seq

    alignment.close()

    return seq_dict

def add_phred_curve(ax, beginning_pos, ending_pos, read_count=20, average=False):
    global reference_seq
    subset_df = pd.DataFrame(index=range(read_count),
                             columns=range(beginning_pos, ending_pos))

    seq = reference_seq

    # Normalize
    for i1, i2 in enumerate(range(read_count)):
        qname = phred_df.iloc[i2, 1]
        quals = phred_df.iloc[i2, 2:].values
        for j1, j2 in enumerate(range(beginning_pos, ending_pos)):
            kmer = get_5_mer(seq, j2).upper()
            avg_phred = avg_phred_per_kmer_df.loc[kmer][0]
            plain_phred = quals[j2]
            subset_df.iloc[i1, j1] = plain_phred - avg_phred

    if average:
        plot_df = subset_df.mean().reset_index(drop=True)
    else:
        plot_df = subset_df.T.reset_index(drop=True)

    plot_df.plot(ax = ax, legend = None)

def plot_phred_curve(beginning_pos, ending_pos, average=False, title="Per Base Sequence Quality", output_file="aligned_plot.png"):
    length = ending_pos - beginning_pos + 1
    fig, ax = plt.subplots(figsize=(12,5))

    # Adding the lines.
    add_curve(ax, phred_df, beginning_pos, ending_pos, average=average)
    add_react_curve(ax, beginning_pos, ending_pos, average=average)
    add_phred_curve(ax, beginning_pos, ending_pos, average=average)

    ax.set_xticks(np.arange(0, length, 5))
    ax.set_xticklabels(np.arange(beginning_pos, ending_pos + 1, 5))
    ax.set_xlabel('position (bp)')
    ax.set_xlim((0, length))
    # ax.set_ylim((-20, 20))

    ax.set_title(title)
    fig.savefig(output_file)
    plt.close(fig)

def create_and_save_phred_df(alignment_file=None, csv_file=None):
    df = read_aligned_sequences(alignment_file)
    df.to_csv(csv_file)

seq_dict = None
phred_df = None
dwell_df = None

###############################################################################
#                                  reactivity                                 #
###############################################################################

react_df = None

def read_pairing_probs(filename, length):
    react_df = pd.DataFrame(index=range(length),
                            columns=["prob"])

    df = pd.read_csv(filename, sep='\t', skiprows=1)

    for row in df.iterrows():
        react_df.loc[row[1].i, "prob"] = row[1]["-log10(Probability)"]

    return react_df

def add_react_curve(ax, beginning_pos, ending_pos, read_count=20, average=False, label=None):
    global react_df
    subset_df = react_df.iloc[beginning_pos:ending_pos]
    plot_df = subset_df.reset_index(drop=True)
    plot_df.plot(ax = ax, label = label)

# react_df = read_pairing_probs("../nmeth.3029-S4.txt", 9173)

###############################################################################
#                             secondary structure                             #
###############################################################################

def parse_nums(s, i, num):
    while s.iloc[i] == num or np.isnan(s.iloc[i]):
        i += 1

        if i >= len(s) - 1:
            break

    return i

def parse_site(s, i):
    i = parse_nums(s, i, 5.0)
    i = parse_nums(s, i, 1.0)
    i = parse_nums(s, i, 3.0)
    return i

def trim_site(s, site):
    i = site[0]
    j = site[1]

    while s.iloc[i] != 5.0:
        i += 1

    while s.iloc[j] != 3.0:
        j -= 1

    return (i, j)

def get_sites(s):
    sites = []
    i = 0

    while i < len(s) - 1:
        j = parse_site(s, i)
        site = trim_site(s, (i, j))
        sites.append(site)
        i = j

    return sites

###############################################################################
#                                  exporting                                  #
###############################################################################

def export(dfs):
    df_names = ["avg_ctrl_phred", "avg_ctrl_dwell", "avg_trizol_dwell", "avg_trizol_phred", "avg_dna_phred", "avg_dna_dwell"]

    col_names = []
    cols = []

    full_excel = "/home/tassos/work/full_structure_annotations.xlsx"
    full_df = pd.read_excel(full_excel)
    sequence = full_df.iloc[1:, 1]
    sequence.index = sequence.index - 1
    sequence.name = "sequence"

    for name, df in zip(df_names, dfs):
        col_names.append(name)
        col = df.iloc[:, 2:].mean()
        cols.append(col)

    col_names.append("secondary structure")
    cols.append(struct)

    df = pd.DataFrame(sequence)
    df.index = cols[0].index

    for name, col in zip(col_names, cols):
        col.name = name
        df = df.join(col)
