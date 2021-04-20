# standard imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
from Bio import SeqIO

# displays logos inline within the notebook;
# remove if using a python interpreter instead
#matplotlib inline
aa_group_list_logo = ['P', 'N', 'R', 'O', 'L']
import logomaker


df = pd.read_csv('data/helix_10_groups.csv')
df = df.sort_values(by=['translated_into_groups'])
df1 = df.iloc[0:200]
df2 = df.iloc[3500:3700]
df = df1.append(df2, ignore_index=True)
num_seqs = len(df)
scores = [[0 for i in range(num_seqs)] for j in range(num_seqs)]


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


for i in range(0, num_seqs):
    for j in range(0, num_seqs):
        score = hamming_distance(df.loc[i, 'translated_into_groups'], df.loc[j, 'translated_into_groups'])
        scores[i][j] = score

num_clusters = 10
print('start kmeans...')
kmeans = cluster.KMeans(num_clusters)
print('ok kmeans...')
results = kmeans.fit(scores)
print('fit kmeans...')
labels = results.labels_
clusters = [[] for i in range(num_clusters)]

for i in range(0,num_seqs):
    clusters[labels[i]].append(df.loc[i, 'translated_into_groups'])

for items in clusters:
    create_logo(items)
    plt.show()
    print('cluster completed')



