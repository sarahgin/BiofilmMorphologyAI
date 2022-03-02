import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import logomaker

from BioinformaticsLab.ComputationalBiology.data_analysis.all_features_calculator import KmerFeatures
from BioinformaticsLab.ComputationalBiology.data_analysis.gene_features_calculator import create_prefix_suffix_dict
from BioinformaticsLab.ComputationalBiology.data_analysis.kmers_analysis import create_prefix_suffix_agg_df


def visualize_prefix_suffix(prefix_suffix_dict, prefix_count):
    # normalize
    for k in prefix_suffix_dict:
        for j in prefix_suffix_dict[k]:
            prefix_suffix_dict[k][j] /= prefix_count[k]
            prefix_suffix_dict[k][j] *= 100
            #if prefix_count[k] < 30:
            #    prefix_suffix_dict[k][j] = 0

    max_dict = {}
    for k in prefix_suffix_dict:
        d = prefix_suffix_dict[k]
        d_max_val = max(d.values())

        # taking list of car values in v
        v = list(d.values())
        # taking list of car keys in v
        vk = list(d.keys())
        d_max_str = (vk[v.index(max(v))])

        max_dict[k + '-' + d_max_str + ' (' + str(prefix_count[k]) + ')'] = d_max_val

    # sorted_max_dict = {k: v for k1, v1 in sorted(max_dict.items(), key=lambda item: item[1])}

    z = list(range(len(max_dict.values())))
    y = list(max_dict.values())

    fig, ax = plt.subplots()
    ax.scatter(z, y)

    for i, txt in enumerate(max_dict.keys()):
        if y[i] > 0:
            ax.annotate(txt, (z[i], y[i]))

    plt.show()


if __name__ == '__main__':
    df = pd.read_csv('data/helix_10_groups.csv')
    df[KmerFeatures.PREFIX_SUFFIX_DICT.name] = df['translated_into_groups'] \
        .apply(lambda helix_seq: create_prefix_suffix_dict(helix_seq, 8, 8, 2, 2))
    agg_dict, count_dict = create_prefix_suffix_agg_df(df)
    visualize_prefix_suffix(agg_dict, count_dict)
    print('done')
