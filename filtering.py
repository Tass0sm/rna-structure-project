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
#                            Average Phred Per Kmer                           #
###############################################################################

bases = ["A", "C", "G", "T"]

k = 5
all_kmers = ["".join(t) for t in product(*[bases] * k)]

def get_5_mer(seq, i):
    return seq[(i - 2):(i + 3)]

def get_phred_scores_for_kmers(filename):
    kmer_dict = dict.fromkeys(all_kmers)

    for kmer in kmer_dict.keys():
        kmer_dict[kmer] = []

    alignment = AlignmentFile(filename)

    for record in tqdm(alignment.fetch()):
        seq = record.seq
        qual = record.qual

        for i in range(2, len(seq) - 3):
            kmer = get_5_mer(seq, i)
            q = ord(qual[i]) - 33
            kmer_dict[kmer].append(q)

    alignment.close()

    return kmer_dict

def get_avg_phred_score_for_kmers(filename):
    kmer_dict = get_phred_scores_for_kmers(filename)
    avg_dict = dict.fromkeys(all_kmers, 0)

    for kmer, qs in tqdm(kmer_dict.items()):
        avg_dict[kmer] = sum(qs) / len(qs)

    return avg_dict

def get_series_5_mer(seq, i):
    j = int(i)
    return seq[(j - 2):(j + 3)].values.sum()

def normalize_by_kmer(curve, sequence, kmer_dict):
    def get_normalized_phred_item(row):
        i = row[0]
        kmer = get_series_5_mer(sequence, i)
        phred = row[1]
        avg_phred = kmer_dict[kmer]
        return phred - avg_phred

    def get_normalized_phred_row(row):
        i = row[0]
        kmer = get_series_5_mer(sequence, i)
        phred_row = row[1:]
        avg_phred = kmer_dict[kmer]
        return phred_row.map(lambda p: p - avg_phred)

    temp_curve = curve.reset_index()

    if isinstance(curve, pd.Series):
        new_curve = pd.Series(temp_curve.apply(get_normalized_phred_item, axis=1).values,
                              index=curve.index)
        return new_curve
    elif isinstance(curve, pd.DataFrame):
        new_curve = pd.DataFrame(temp_curve.apply(get_normalized_phred_row, axis=1, result_type="expand").values,
                                 index=curve.index)
        return new_curve
    else:
        return None
