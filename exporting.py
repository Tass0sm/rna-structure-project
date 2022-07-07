import pandas as pd

default_df_names = ["avg_ctrl_phred", "avg_ctrl_dwell", "avg_trizol_dwell", "avg_trizol_phred", "avg_dna_phred", "avg_dna_dwell"]

def export(dfs, df_names):
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
