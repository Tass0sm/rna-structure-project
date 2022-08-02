import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess

from Bio import SeqIO

from itertools import product
from pysam import AlignmentFile
from tqdm import tqdm

import utils

###############################################################################
#                             Getting Phred Scores                            #
###############################################################################

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
#                                     Old                                     #
###############################################################################

def align_fastq(filename, output):
    reference = "/home/tassos/work/common-data/RNA_section__454_9627.fasta"
    script = "/home/tassos/work/rna-structure-project/scripts/align_reads.sh"
    subprocess.run([script, reference, filename, output])

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
