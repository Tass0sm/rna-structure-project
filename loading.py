import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess

import csv
from tqdm import tqdm
import numpy as np

from Bio import SeqIO

from itertools import product
from pysam import AlignmentFile
from tqdm import tqdm

import utils

# def read_seq_fasta_and_combine_with_average_phred():
#     dna_seq_df = pd.read_csv("./data/dna_seq_df.csv")
#     dna_seq_df = dna_seq_df.set_index("Unnamed: 0")
#     mean_phred_df = pd.DataFrame(reverse_dna_phred_df.mean().iloc[1:]).T
#     combined_df = pd.concat([dna_seq_df, mean_phred_df], axis=0)

# def read_alignment_and_save():
#     o = f"./data/dna_sense/alignment.sam"
#     df = loading.read_aligned_sequences(o)
#     c = f"./data/dna_sense/dna_sense_phreds.csv"
#     df.to_csv(c)

###############################################################################
#                             Getting Phred Scores                            #
###############################################################################

def get_rna_seq_df(filename):
    result = dict()

    for record in SeqIO.parse("/fs/project/PAS1405/General/HIV_RNA_modification_dataset/RNA_section__454_9627.fasta", "fasta"):
        seq = record.seq

    s = pd.Series(iter(seq))
    df = pd.DataFrame(s).T
    return df

def get_reference_seq_dict():
    result = dict()

    for record in SeqIO.parse("./data/human/gencode.v36.transcripts.fa", "fasta"):
        result[record.id] = record.seq

    return result

def read_alignment_with_multiple_references():
    ref_dict = get_reference_seq_dict()

    # read_len = alignment.count()
    # alignment.close()

    # col_index = ["read_id"] + list(range(ref_len))

    # results_df = pd.DataFrame(index=range(read_len),
    #                           columns=col_index)

    alignment_file = AlignmentFile("./data/human/alignment.sam")
    reference_name_to_phreds_dict = dict()

    for i, record in tqdm(enumerate(alignment_file.fetch())):
        qual_seq = record.query_qualities

        if qual_seq is not None:
            ref_name = record.reference_name
            ref_seq = ref_dict[ref_name]
            ref_len = len(ref_seq)

            read_seq = record.query_sequence
            alignment = record.get_reference_positions(full_length=True)

            result_array = np.full(ref_len, np.nan)

            for p, q in zip(alignment, qual_seq):
                if (p is not None):
                    result_array[p] = q

            result_series = pd.Series(result_array)
            result_series.name = record.query_name

            if ref_name not in reference_name_to_phreds_dict:
                reference_name_to_phreds_dict[ref_name] = []

            reference_name_to_phreds_dict[ref_name].append(result_series)

    alignment_file.close()

    return reference_name_to_phreds_dict

rs = "ENST00000224237.9|ENSG00000026025.16|OTTHUMG00000017744|OTTHUMT00000047015.1|VIM-201|VIM|1868|protein_coding|"
best_r = 'ENST00000387347.2|ENSG00000210082.2|-|-|MT-RNR2-201|MT-RNR2|1559|Mt_rRNA|'
# https://genome-asia.ucsc.edu/cgi-bin/hgc?hgsid=771583428_tT7LntFBfwW5yoyaJ42HibfdgbLR&g=htcDnaNearGene&i=ENST00000387347.2&c=chrM&l=1670&r=3229&o=knownGene&boolshad.hgSeq.promoter=0&hgSeq.promoterSize=1000&hgSeq.utrExon5=on&boolshad.hgSeq.utrExon5=0&hgSeq.cdsExon=on&boolshad.hgSeq.cdsExon=0&hgSeq.utrExon3=on&boolshad.hgSeq.utrExon3=0&hgSeq.intron=on&boolshad.hgSeq.intron=0&boolshad.hgSeq.downstream=0&hgSeq.downstreamSize=1000&hgSeq.granularity=gene&hgSeq.padding5=0&hgSeq.padding3=0&boolshad.hgSeq.splitCDSUTR=0&hgSeq.casing=exon&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&submit=submit

def read_phreds_for_reference(ref_name):
    ref_dict = get_reference_seq_dict()
    ref_seq = ref_dict[ref_name]
    ref_len = len(ref_seq)

    col_index = ["read_id"] + list(range(ref_len))
    result_series_list = []

    alignment_file = AlignmentFile("./data/human/alignment.sam")

    for i, record in tqdm(enumerate(alignment_file.fetch())):
        # print(record.reference_name)
        if record.reference_name == ref_name:
            qual_seq = record.query_qualities

            if qual_seq is not None:
                read_seq = record.query_sequence
                alignment = record.get_reference_positions(full_length=True)

                result_array = np.full(ref_len, np.nan)

                for p, q in zip(alignment, qual_seq):
                    if (p is not None):
                        result_array[p] = q

                result_series = pd.Series(result_array)
                result_series.name = record.query_name
                result_series_list.append(result_series)

    alignment_file.close()

    result_df = pd.concat(result_series_list, axis=1).T

    return result_df

def count_reads_for_each_reference():
    alignment_file = AlignmentFile("./data/human/alignment.sam")
    reference_name_to_count_dict = dict()

    for i, record in tqdm(enumerate(alignment_file.fetch())):
        ref_name = record.reference_name

        if ref_name not in reference_name_to_count_dict:
            reference_name_to_count_dict[ref_name] = 0

        reference_name_to_count_dict[ref_name] += 1

    alignment_file.close()

    return reference_name_to_count_dict

def highest_read_references():
    d = count_reads_for_each_reference()
    d_sorted_list = sorted(d.items(), key=lambda x:x[1], reverse=True)
    new_d = dict(d_sorted_list)
    return new_d


###############################################################################
#                                     csvs                                    #
###############################################################################

def load_and_avg_csv(filename):
    col_len = 0

    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')

        for r, row in tqdm(enumerate(reader)):
            if r == 0:
                col_len = len(row)
                print(f"Len: {col_len}")
                row_array = np.zeros(col_len)
                count_array = np.zeros(col_len)

            for c in range(0, col_len - 1):
                if row[c + 1] != "":
                    row_array[c] += float(row[c + 1])
                    count_array[c] += 1

    return np.divide(row_array, count_array)

###############################################################################
#                                     Old                                     #
###############################################################################

def align_fastq(filename, output):
    reference = "/users/PAS1405/tassosm/Desktop/common-data/RNA_section__454_9627.fasta"
    script = "./scripts/align_reads.sh"
    subprocess.run([script, reference, filename, output])

def read_aligned_sequences(filename, limit=None):
    alignment = AlignmentFile(filename)
    ref_len = alignment.lengths[0]
    read_len = alignment.count()
    alignment.close()

    col_index = ["read_id"] + list(range(ref_len))

    if limit is not None:
        read_len = limit + 2

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

        if limit is not None and i > limit:
            break

    alignment.close()

    return results_df

def read_and_average_aligned_sequences(filename):
    alignment = AlignmentFile(filename)
    ref_len = alignment.lengths[0]
    read_len = alignment.count()
    alignment.close()

    alignment = AlignmentFile(filename)

    row_array = np.zeros(ref_len)
    count_array = np.zeros(ref_len)

    for i, record in tqdm(enumerate(alignment.fetch())):
        if (record.qual is not None):
            for p, q in zip(record.get_reference_positions(full_length=True),
                            record.qual):
                if (p is not None):
                    phred = ord(q) - 33
                    row_array[p] += phred
                    count_array[p] += 1

    alignment.close()

    avg_array = np.divide(row_array, count_array)
    return avg_array
