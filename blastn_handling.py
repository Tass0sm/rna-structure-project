import os
import pandas as pd

def filter_and_output_subread_lists():
    base_dir = "/fs/ess/PAS1405/ForTassos/"
    dirs = ["01012022_LambdaFragment_RCA", "01152021_LambdaFragment_RCA", "01152021_LambdaHindIII_RCA", "05182021_LambdaHindIII_RCA", "12072021+12142021_Ecoli_RCA"]
    seq_csvs = [os.path.join(base_dir, d, "df_w_seqs.csv") for d in dirs]
    subread_csvs = [os.path.join(base_dir, d, "subread_blastn_ws7.csv") for d in dirs]

    def get_direction_proportion(g):
        return (g.direction == "plus").mean()

    for d_name, seq_csv, subread_csv in zip(dirs, seq_csvs, subread_csvs):
        seqs_df = pd.read_csv(seq_csv)
        subread_df = pd.read_csv(subread_csv, sep="\t")
        subread_df.columns = ["name", "direction"]

        seqs_df_with_direction = seqs_df.merge(subread_df, on="name")
        seqs_df_with_direction.drop(columns=["Unnamed: 0"], inplace=True)
        grouped_seqs_df = seqs_df_with_direction.groupby(by="qseqid")

        # filter by number of subreads
        grouped_seqs_df = grouped_seqs_df.filter(lambda g: g.shape[0] > 10)
        grouped_seqs_df = grouped_seqs_df.groupby(by="qseqid")

        # filter by direction ratio
        grouped_seqs_df = grouped_seqs_df.filter(lambda g: 0.4 <= get_direction_proportion(g) <= 0.6)
        grouped_seqs_df = grouped_seqs_df.groupby(by="qseqid")

        # output_df
        output_df = grouped_seqs_df.apply(lambda g: pd.Series([g.shape[0], get_direction_proportion(g)]))
        output_df.columns = ["number_of_subreads", "proportion_of_plus_subreads"]
        output_df.to_csv(f"./{d_name}_output.csv")

