import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

from Bio import Phylo
from matplotlib_venn import venn2, venn3
import itertools
import scipy
from scipy.spatial import distance
import pickle

from ComputationalBiology.biore.biore_utils import aa_into_group, aa_into_BY, get_hamming
from ComputationalBiology.biore.visualize_utils import create_logo

pd.set_option('display.max_columns', None)


def parse_file(file_path):
    """
    The function parses the given tab-separated file and returns a dataframe object
    :return:
    """
    df = pd.read_csv(file_path, sep='\t', header=0)
    return df


def get_subsequences(motif_info: str, sequence: str, prefix_type: str):
    """
    Find all subsequences in sequence according to the ranges defines in motif_info
    :param motif_info: String describing motif(s) locations within sequence
    :param sequence: query sequence
    :return: a list of subsequences
    """
    res = []
    # find all matches of (HELIX (\d+)..(\d+))# TODO: make general fo beta
    all_matches = re.findall(prefix_type + ' (\d+)\\.\\.(\d+)', motif_info)
    assert (all_matches is not None)

    # go over the matches
    for start, end in all_matches:
        # print(int(start)- int(end))
        res.append(sequence[int(start) - 1:int(end) + 1 - 1])

    return res


def create_subseq_df(df: pd.DataFrame, col_name: str, subseq_type: str):
    """
    Given a dataframe where each row is a gene (protein) we output a new dataframe
    where each row is a subsequence (e.g. helix, transmembrane etc) with another
    column of the matching length of the subsequence.
    :param df: Input df
    :param col_name: Name of the column in the original df (e.g. 'Helix', 'Transmembrane')
    :param subseq_type: e.g. 'HELIX' \ 'TRANSMEM' etc.
    :return:
    """

    df2 = df[~df[col_name].isna()].copy()  # remove all rows which do not have the subseq

    # parse the content of {col_name} to fetch a list of all subsequences
    df2['seq_' + col_name] = \
        df2.apply(lambda row: get_subsequences(row[col_name], row['Sequence'], subseq_type), axis=1)

    df3 = df2[['seq_' + col_name, 'Protein names']].copy()

    print('num of genes with at least one subseq: ', len(df3))
    # split each row into multiple rows - according to the number of found subsequences
    df_all_sequences = pd.DataFrame(df3.explode(column='seq_' + col_name))

    # add the length row
    df_all_sequences['subseq_len'] = df_all_sequences['seq_' + col_name].apply(lambda x: len(x))
    # df_all_sequences['subseq_len'].plot.hist(bins=100)
    print('num of total subseqs: ', len(df_all_sequences))
    return df_all_sequences


def find_S_by_column(df: pd.DataFrame, col_name: str):
    df_counter = df.groupby(col_name).count()
    df_counter['sequence'] = df_counter['Protein names'] #since sequence column is now index
    assert (len(df[col_name].unique() == len(df_counter)))

    max_count = df_counter['sequence'].max()
    max_value = df_counter[df_counter['sequence'] == max_count].index.tolist()[0]
    total_seqs_non_unique = sum(df_counter['sequence'])
    print('col_name={}: {} was found: {} ({:.2}%),out of {} valid seqs'
          .format(col_name, max_value, max_count, 100 * max_count / total_seqs_non_unique, total_seqs_non_unique))
    # df_counter['sequence'].plot.hist(bins=2500)
    assert sum(df[col_name] == max_value) == max_count  # sainty check
    return df_counter


def plot_cumsum(df_counter, color, figHandle, columm_name, to_show=False):
    max_count = df_counter[columm_name].nlargest(len(df_counter))
    percent_values = max_count.values / sum(max_count) * 100
    y_cumsum = np.cumsum(percent_values)

    x_cumsum = np.arange(1, len(df_counter)+1)/len(df_counter)

    plt.plot(x_cumsum, y_cumsum, 'o', figure=figHandle)
    plt.xlabel('Percentage of sequences')
    plt.ylabel('Cumulative percentage covered by sequences')
    if to_show:
        plt.show()

def get_all_substrings(sequence, min_len):
    res = [sequence[i: j] for i in range(len(sequence)) for j in range(i + 1, len(sequence) + 1) if
           len(sequence[i:j]) == min_len]
    return res


def generate_venns():
    humanFile = open('data//human1//top_subseqs_human1.txt', 'r')
    humanSet = humanFile.readlines()
    humanSetUnique = np.unique(humanSet)

    mouseFile = open('data//mouse1//top_subseqs_mouse1.txt', 'r')
    mouseSet = mouseFile.readlines()
    mouseSetUnique = np.unique(mouseSet)

    intersectionSet = set(humanSetUnique) & set(mouseSetUnique)

    print(len(humanSet), ',', len(mouseSet), ',')
    print(len(humanSetUnique), ',', len(mouseSetUnique), ',', len(intersectionSet))

    venn2(subsets=(len(humanSetUnique), len(mouseSetUnique), len(intersectionSet)),
          set_labels=('Human', 'Mouse'))
    plt.show()


def organisms_analysis(organism: str):
    dfs = []

    file_path = 'data/' + organism + '.tab'
    df = parse_file(file_path)

    #split into reviewed vs. unreviewed
    #leave 'reviewed' only
    reviewedDF = df[df['Status'] == 'reviewed']
    unreviewedDF = df[df['Status'] == 'unreviewed']
    assert (len(df) == len(reviewedDF) + len(unreviewedDF))
    df = reviewedDF
    print('num of genes in raw data: ', len(df),
          ' Reviewed: ', len(reviewedDF),
          'Unreviewed: ', len(unreviewedDF))

    #create df with onlu transmembrane
    df_transmembrane = create_subseq_df(df=df, col_name='Transmembrane', subseq_type='TRANSMEM').copy()
    df_transmembrane = df_transmembrane.rename(columns={"seq_Transmembrane": "sequence"})

    if False:  # For each length, the total number of subsequences at that length
        df_transmembrane['subseq_len'].plot.hist(bins=np.arange(1, 50),
                                                 alpha=0.5,
                                                 figure=fig,
                                                 label=organism)
        plt.legend()

    min_len = 21
    max_len = 22

    for j, chosen_subseq_len in enumerate(range(min_len, max_len)):
        #if chosen_subseq_len != 21:
        #    continue

        print('Analyzing len: ', str(chosen_subseq_len))
        df_curr_length = df_transmembrane[df_transmembrane['subseq_len'] == chosen_subseq_len].copy()
        print('num of subseqs of length {}: {} '.format(chosen_subseq_len, len(df_curr_length)))

        # remove sequences with invalid amino acids
        df_curr_length['is_valid_seq'] = df_curr_length['sequence'].apply(
            lambda seq: 'X' not in seq).copy()

        df_curr_length = df_curr_length[df_curr_length['is_valid_seq']].copy()
        n_total_unique = len(df_curr_length)
        if n_total_unique == 0:
            continue
        print('num of valid subseqs of length {}: {} '.format(chosen_subseq_len, len(df_curr_length)))
        print('num unique subseqs:', len(df_curr_length['sequence'].unique()))

        #organize and sort by unique sequences
        df_counter_all = find_S_by_column(df_curr_length, 'sequence')
        df_counter_all_sorted = df_counter_all.sort_values(by='subseq_len', ascending=False)

        dfs.append(df_counter_all_sorted)
        print('----Finished current length: ', str(chosen_subseq_len), '-----')

    print('done all lengths')
    return dfs


def parse_IEDB_excel():
    filenames = ['data//human1//human1.0.csv',
                 'data//human1//human1.1.csv',
                 'data//human1//human1.2.csv',
                 'data//human1//human1.3.csv',
                 'data//human1//human1.4.csv']
    total_consensus = 0
    merged_dict = {}  # will contain consensus as keys and matching values as the list of sequences
    singletons_merged = []
    for filename in filenames:
        print(filename + '...')
        df = pd.read_csv(filename, dtype=str)

        df_non_consensus = df[df['Peptide Number'] != 'Consensus']
        print('Total sequences: ', len(df_non_consensus))

        df_singletons = df[df['Peptide Number'] == 'Singleton']
        singletons_merged.append(df_singletons)
        print('Total singletons: ', len(df_singletons))

        df_concensus = df[df['Peptide Number'] == 'Consensus']
        print('Total consensus sequences: ', len(df_concensus))
        total_consensus += len(df_concensus)

        #only sequences remain
        df = df[(df['Peptide Number'] != 'Consensus') &
                (df['Peptide Number'] != 'Singleton')]

        #for each cluster get the consensus string and write
        for cluster_id in set(df['Cluster.Sub-Cluster Number'].values):
            current_df = df[df['Cluster.Sub-Cluster Number'] == cluster_id].copy()
            cluster_sequences = current_df['Peptide'].values
            cluster_consensus = df_concensus[df_concensus['Cluster.Sub-Cluster Number'] == cluster_id]['Alignment'].values

            # adding a new entry to merged_dict (all 4 files)
            assert len(cluster_consensus) == 1
            cluster_consensus_seq = cluster_consensus[0]
            if cluster_consensus_seq not in merged_dict.keys():
                merged_dict[cluster_consensus_seq] = list(cluster_sequences)
            else: # i.e. this consensus exists already
                merged_dict[cluster_consensus_seq].append(list(cluster_sequences))

    print('All files: len of df_merged: {} vs. len of df_concensus: {}'.
          format(len(merged_dict), total_consensus))

    #create df_merged: concensus_str, [list of subseqs], #subseqs, #X in concensus
    df_merged = create_df_from_dict(merged_dict)

    return df_merged, pd.concat(singletons_merged)


def create_df_from_dict(data_dict):
    # data_dict: key is sequence, value is the list of covered sequence(s)
    data = []
    for concensus_str in data_dict:
        data.append({'concensus': concensus_str, 'subseqs':data_dict[concensus_str]})
    df_res = pd.DataFrame(data)
    df_res['num_seqs'] = df_res['subseqs'].apply(lambda x: len(x))
    df_res['num_X'] = df_res['concensus'].apply(lambda x: x.count('X'))
    df_res['len_concensus'] = df_res['concensus'].apply(lambda x: len(x))
    return df_res


def add_singeltons_to_df_merged(df_merged, df_singletons_merged):
    # create a temp df from df_singletons_merged in order to concatenate the two dfs

    keys = df_singletons_merged['Alignment'].values
    # values = [[k] for k in keys]
    # singletons_dict = dict(zip(keys, values))

    singletons_dict = dict((key, [key]) for key in keys)

    df_temp = create_df_from_dict(singletons_dict)
    return pd.concat([df_merged, df_temp])


def fill_missing_seqs(df_human, df_human_merged):
    # df_human is the original data where each unique sequence could appear more tan once,
    # in df_human_merged each sequence apeared once only
    # Goal: fill the remaining once such that the total covered sequences will be sum(df_human['subseq_len'])

    # for each sequence in df_human that appears n > 1 times:
    # find its cluster
    # add it n-1 times
    # TODO
    return 0


def create_lower_triang_distance_matrix(sequences: list):
    dist_matrix = np.zeros((len(sequences), len(sequences)))
    for i in range(0, len(sequences)):
        c1 = sequences[i]
        if i % 50 == 0:
            print(i)

        for j in range(i+1, len(sequences)):
            c2 = sequences[j]

            # TODO: use a "smarter" distance matrix for peptide comparison (e.g. our 5-groups or BLOSUM62)
            dist_matrix[j][i] = get_hamming(c1, c2)

    return dist_matrix


def create_triangle_distance_matrix(matrix):
    res = []
    temp = np.tril(dist_matrix).tolist()

    for i in range(len(temp)):
        res.append(temp[i][0:i+1])

    return res



if __name__ == '__main__':

    #STEP 1 - generate all sequences of length 21
    dfs = organisms_analysis('human1')
    df_human = dfs[0]
    for l in range(len(dfs)):
        fig = plt.figure()
        plot_cumsum(dfs[l], 'r', fig, columm_name='sequence', to_show=False)

    #STEP 2 - split 15K into 5x3K
    #STEP 3 - run IEDB 5 times (http://tools.iedb.org/cluster/)
    #STEP 4 - parse IEDB results csv files
    df_merged, df_singletons_merged = parse_IEDB_excel()
    # TODO: fill_missing_seqs

    # STEP 5 - add singletons to the df
    df_merged_with_singletons = add_singeltons_to_df_merged(df_merged, df_singletons_merged)
    plot_cumsum(df_merged_with_singletons, 'r', fig, columm_name='num_seqs', to_show=False)
    plt.show()
    # assert (sum(d1['subseq_len']) == sum(df_merged_with_singletons['num_seqs'])) % assertion after running fill_missing_seqs

    # Create distance matrix for all 15333 sequences

    # dist_matrix = create_lower_triang_distance_matrix(list(df_human.index.values))
    # all_concensuses =list(df_merged[df_merged['len_concensus'] == 21]['concensus'].values)
    all_sequences = list(df_human.index.values)
    dist_matrix = create_lower_triang_distance_matrix(all_sequences)
    dist_lists = create_triangle_distance_matrix(dist_matrix)

    M_file = './dist_lists.pickle'
    with open(M_file, 'wb') as pickle_file:
        pickle.dump(dist_lists, file=pickle_file)

    with open(M_file, 'rb') as pickle_file:
        dist_lists = pickle.load(file=pickle_file)

    # Create a tree
    import Bio
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
    constructor = DistanceTreeConstructor()
    print('Constructing tree using UPGMA')
    tree = constructor.upgma(DistanceMatrix(names=all_sequences, matrix=dist_lists))

    tree_file = './tree.pickle'
    with open(tree_file, 'wb') as pickle_file:
        pickle.dump(tree, file=pickle_file)

    # with open(tree_file, 'rb') as pickle_file:
    #     tree2 = pickle.load(file=pickle_file)

    Bio.Phylo.draw(tree)
    print('done')

    # M2[M2 >= 0.5] = 1
    # M2[M2 < 0.5] = 0
    # import seaborn as sns
    # sns.heatmap(M2)
    # plt.show()


    # assert(np.all(np.equal(M, M2)))
    print('done')


    # generate_venns()
