# standard imports
import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
import logomaker


def hamming_distance(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


def create_logo(items):
    df = pd.DataFrame(items, columns=['translated_into_groups'])

    for p in range(10):
        df['position_' + str(p)] = df['translated_into_groups'].apply(lambda x: x[p])

    df_logo = pd.DataFrame(np.zeros((10, 5)), columns=aa_group_list_logo)

    for p in range(10):
        for g in aa_group_list_logo:
            df_logo.loc[p, g] = (sum(df['position_' + str(p)] == g))

    seq_num = df_logo.sum(axis=1).iloc[0]
    df_logo = df_logo.applymap(lambda i: i / seq_num)

    ss_df = df_logo
    ss_logo = logomaker.Logo(ss_df,
                             width=.8,
                             vpad=.05,
                             fade_probabilities=True,
                             stack_order='small_on_top',
                             color_scheme='dodgerblue',
                             font_name='Arial')

    # style using Logo methods
    ss_logo.style_spines(spines=['left', 'right'], visible=False)

    # style using Axes methods
    ss_logo.ax.set_xticks(range(len(ss_df)))
    ss_logo.ax.set_xticklabels('%+d' % x for x in [-3, -2, -1, 1, 2, 3, 4, 5, 6])
    ss_logo.ax.set_yticks([0, .5, 1])
    ss_logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')
    ss_logo.ax.set_ylabel('probability')
    plt.show()


# this main does kmeans over 10-length helix sequences (that were translated into groups)
if __name__ == '__main__':
    aa_group_list_logo = ['P', 'N', 'R', 'O', 'L']
    df = pd.read_csv('data/helix_10_groups.csv')
    df = df.sort_values(by=['translated_into_groups'])

    #df1 = df.iloc[0:100]
    #df2 = df.iloc[3500:3600]
    #df = df1.append(df2, ignore_index=True)

    num_seqs = len(df)

    SCORES_PICKLE_FILE = './data/helix_hamming_scores.pickle'
    if not os.path.exists(SCORES_PICKLE_FILE):
        # compute scores
        scores = [[0 for i in range(num_seqs)] for j in range(num_seqs)]
        for i in range(0, num_seqs):

            if i % 100 == 0:  # printing for debug only
                print(i, ' out of ', num_seqs)

            for j in range(0, num_seqs):
                score = hamming_distance(df.loc[i, 'translated_into_groups'], df.loc[j, 'translated_into_groups'])
                scores[i][j] = score
        # save to pickle
        with open(SCORES_PICKLE_FILE, 'wb') as pickle_file:
            pickle.dump(scores, file=pickle_file)
    else:
        with open(SCORES_PICKLE_FILE, 'rb') as pickle_file:
            scores = pickle.load(file=pickle_file)
    if True:
        # k-means clustering of helix sequences
        num_clusters = 10
        kmeans = cluster.KMeans(num_clusters)
        results = kmeans.fit(scores)
        labels = results.labels_
        clusters = [[] for i in range(num_clusters)]
        for i in range(0, num_seqs):
            clusters[labels[i]].append(df.loc[i, 'translated_into_groups'])

        # create logo for each cluster
        for items in clusters:
            create_logo(items)
            plt.show()
