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

import utils

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
