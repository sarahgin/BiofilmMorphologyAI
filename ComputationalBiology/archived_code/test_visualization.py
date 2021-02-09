import matplotlib.pyplot as plt
import pandas as pd

from ComputationalBiology.data_analysis.main_parser_features_calc import NEXT_NT_DF_FILE

if __name__ == '__main__':
    print('Loading pickle file: {}...'.format(NEXT_NT_DF_FILE))
    next_nt_df = pd.read_pickle(NEXT_NT_DF_FILE)
    next_nt_df["SUM"] = next_nt_df[["A", "C", "G", "T"]].sum(axis=1)
    next_nt_df["A_PER"] = next_nt_df["A"] / next_nt_df["SUM"] * 100
    next_nt_df["C_PER"] = next_nt_df["C"] / next_nt_df["SUM"] * 100
    next_nt_df["G_PER"] = next_nt_df["G"] / next_nt_df["SUM"] * 100
    next_nt_df["T_PER"] = next_nt_df["T"] / next_nt_df["SUM"] * 100

    next_nt_df["MAX_PER"] = next_nt_df[["A_PER", "C_PER", "G_PER", "T_PER"]].max(axis=1)
    plt.scatter(range(len(next_nt_df["MAX_PER"])), next_nt_df["MAX_PER"])
    plt.show()
    exit(0)
