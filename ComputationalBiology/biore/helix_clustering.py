# standard imports
import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import logomaker

helix_length = 10

cluster_by_column = 'translated_into_groups'
aa_group_list_logo = ['B', 'Y']

#cluster_by_column = 'translated_into_groups'
#aa_group_list_logo = ['P', 'N', 'R', 'O', 'L']

#cluster_by_column = 'sequences_Helix'
#aa_group_list_logo = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
alphabet_size = len(aa_group_list_logo)


def hamming_distance(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


def create_logo(items):
    df = pd.DataFrame(items, columns=[cluster_by_column])

    for p in range(helix_length):
        df['position_' + str(p)] = df[cluster_by_column].apply(lambda x: x[p])

    df_logo = pd.DataFrame(np.zeros((helix_length, alphabet_size)), columns=aa_group_list_logo)

    for p in range(helix_length):
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
    return ss_logo


def search_num_cluster(scores):
    """
    :param scores: Data for running K-means, in this case a list of lists
    """
    # run kmeans for different values of K:

    errors = []  # vector for errors
    k_range = range(2, 200)
    for k in k_range:
        print(k)
        kmeans = KMeans(n_clusters=k)
        # kmeans = KMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
        kmeans.fit(scores)  # kmeans algorithm fits to the X dataset

        errors.append(kmeans.inertia_)

        # label = kmeans.labels_
        # from sklearn.metrics import silhouette_score
        # sil_coeff = silhouette_score(scores, label, metric='euclidean')

        # kmeans inertia_ attribute is:  Sum of squared distances of samples to their closest cluster center.

    # Plot the elbow graph
    plt.plot(k_range, errors)
    plt.title('The Elbow Method Graph')
    plt.xlabel('Number of clusters')
    plt.ylabel('errors')
    plt.show()


# this main does kmeans over 10-length helix sequences (that were translated into groups)
if __name__ == '__main__':
    df = pd.read_csv('data/helix_10_groups.csv')
    df = df.sort_values(by=[cluster_by_column])

    # df1 = df.iloc[0:100]
    # df2 = df.iloc[3500:3600]
    # df = df1.append(df2, ignore_index=True)

    num_seqs = len(df)

    SCORES_PICKLE_FILE = './data/helix_hamming_scores.pickle'
    if not os.path.exists(SCORES_PICKLE_FILE):
        # compute scores
        scores = [[0 for i in range(num_seqs)] for j in range(num_seqs)]
        for i in range(0, num_seqs):

            if i % 100 == 0:  # printing for debug only
                print(i, ' out of ', num_seqs)

            for j in range(0, num_seqs):
                score = hamming_distance(df.loc[i, cluster_by_column], df.loc[j, cluster_by_column])
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
        kmeans = KMeans(num_clusters, init='k-means++')
        results = kmeans.fit(scores)
        labels = results.labels_
        clusters = [[] for i in range(num_clusters)]
        for i in range(0, num_seqs):
            clusters[labels[i]].append(df.loc[i, cluster_by_column])

        # create logo for each cluster
        c = 1
        for items in clusters:
            logo = create_logo(items)
            plt.title('Number of helices in this cluster: ' + str(len(items)))
            logo.fig.savefig('./data/cluster_' + str(c) + '.png')
            c = c + 1
