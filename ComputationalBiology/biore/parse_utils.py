# goal: get human1.tab -->  fetch all helix
# save in dataframe
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib_venn import venn2, venn3

from ComputationalBiology.biore.biore_utils import aa_into_group, aa_into_BY


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
    assert (len(df[col_name].unique() == len(df_counter)))

    max_count = df_counter['sequence'].max()
    max_value = df_counter[df_counter['sequence'] == max_count].index.tolist()[0]
    total_seqs_non_unique = sum(df_counter['sequence'])
    print('col_name={}: {} was found: {} ({:.2}%),out of {} valid seqs'
          .format(col_name, max_value, max_count, 100 * max_count / total_seqs_non_unique, total_seqs_non_unique))
    # df_counter['sequence'].plot.hist(bins=2500)
    assert sum(df_transmem[col_name] == max_value) == max_count  # sainty check
    return df_counter


def plot_cumsum(df_counter, color, figHandle):
    max_count = df_counter['sequence'].nlargest(len(df_counter))
    percent_values = max_count.values / sum(max_count) * 100
    y_cumsum = np.cumsum(percent_values)
    plt.plot(y_cumsum, 'o', color=(color, 0, 0), figure=figHandle)


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


if __name__ == '__main__':
    #generate_venns()
    # parsing file to create all input data:
    organisms = ['human1', 'mouse1', 'bsubtilis']
    organismsColors = ['r', 'b', 'g']
    sharedFig = plt.figure()
    cusumFig = plt.figure()
    for o, organism in enumerate(organisms):
        print('Organism: ', organism)
        file_path = 'data/' + organism + '.tab'
        df = parse_file(file_path)

        print('num of genes in raw data: ', len(df))

        df_transmembrane = create_subseq_df(df=df, col_name='Transmembrane', subseq_type='TRANSMEM').copy()
        df_transmembrane['sequence'] = df_transmembrane['seq_Transmembrane']

        if False:  # For each length, the total number of subsequences at that length
            df_transmembrane['subseq_len'].plot.hist(bins=np.arange(1, 50),
                                                     alpha=0.5,
                                                     figure=sharedFig,
                                                     label=organism)
            plt.legend()

        legend_str = []
        min_len = 21
        max_len = 22
        num_of_sequences = 2000
        colors = np.linspace(0, 1, max_len - min_len)

        out_file_aa = open('data//' + organism + '//top_subseqs_' + organism + '.txt', 'w')

        for j, chosen_subseq_len in enumerate(range(min_len, max_len)):
            if chosen_subseq_len != 21:
                continue

            print('Analyzing len: ', str(chosen_subseq_len))
            df_transmem = df_transmembrane[df_transmembrane['subseq_len'] == chosen_subseq_len].copy()
            # Note: there are about 20K subseqs of length: 19, 20, 22, 23 and 36K of 21 bases.

            print('num of subseqs of length {}: {} '.format(chosen_subseq_len, len(df_transmem)))

            # remove sequences with invalid amino acids
            df_transmem['is_valid_seq'] = df_transmem['seq_Transmembrane'].apply(
                lambda seq: 'X' not in seq).copy()

            df_transmem = df_transmem[df_transmem['is_valid_seq']].copy()
            n_total_unique = len(df_transmem)
            if n_total_unique == 0:
                continue
            print('num of valid subseqs of length {}: {} '.format(chosen_subseq_len, len(df_transmem)))

            df_transmem['translated_into_groups'] = df_transmem['seq_Transmembrane'].apply(lambda x: aa_into_group(x))
            df_transmem['translated_into_BY'] = df_transmem['seq_Transmembrane'].apply(lambda x: aa_into_BY(x))

            df_transmem = df_transmem.sort_values(by='translated_into_BY')
            df_transmem['translated_into_BY'].unique()

            print('num unique subseqs:', len(df_transmem['seq_Transmembrane'].unique()))
            #print('num unique 5-groups:', len(df_transmem['translated_into_groups'].unique()))
            #print('num unique BY:', len(df_transmem['translated_into_BY'].unique()))

            df_counter_all = find_S_by_column(df_transmem, 'seq_Transmembrane')

            df_counter_all_sorted = df_counter_all.sort_values(by='subseq_len', ascending=False)

            if False:  # bar plot to show top 1 sequence coverage (the number it covers)
                barX = "[" + str(chosen_subseq_len) + "] "
                barY = df_counter_all_sorted['subseq_len'].head(1).values[0]
                plt.bar(barX, barY * 100 / n_total_unique, figure=sharedFig, color=organismsColors[o], alpha=0.5)

            if False:  # elbow plot
                plot_cumsum(df_counter_all_sorted, colors[j], cusumFig)
                plt.title(organism)

            for TM_sequence_index, s in enumerate(df_counter_all_sorted['subseq_len'].head(num_of_sequences).index):
                all_substrings = get_all_substrings(s, min_len)
                # assert(len(all_substrings) == 1)
                for TM_subsequence_index, subs in enumerate(all_substrings):
                    # out_file_aa.write(">" + organism + "_" + subs + "_"
                    #                  + str(TM_sequence_index)
                    #                  + "[" + str(TM_subsequence_index) + "]"
                    #                  + "\n")           # fasta header
                    out_file_aa.write(subs + "\n")

            top_X_coverage = df_counter_all_sorted['subseq_len'].head(num_of_sequences).sum() / n_total_unique * 100
            print('top ', str(num_of_sequences), ' coverage: ', top_X_coverage)

            # print the top 10- most abundant sequences of current length
            max_count = df_counter_all['sequence'].nlargest(num_of_sequences)
            # print(max_count.index.tolist())

            print('----Finished current length: ', str(chosen_subseq_len), '-----')
            legend_str.append(str(chosen_subseq_len))  # + str(len(df_transmem)))

        out_file_aa.close()

    plt.show()
    print('done')

# if __name__ == '__main__':
#
#     if False:
#         file_path = 'data/human1.tab'
#         df = parse_file(file_path)
#
#         df_all_helices = motif_sequences_to_file(df=df, col_name='Helix', type='HELIX')
#
#         df_all_helices.to_csv('data/all_helix.csv', header=True, sep='\t', index=False)
#
#     #read helix data:
#     print('reading...')
#     df2 = pd.read_csv('data/all_helix.csv')
#     df2['helix_len'] = df2['sequences_Helix'].apply(lambda x: len(x))
#     df2.describe()
#
#     # get only 10-mers
#     df3 = df2.copy()
#     df3 = df3[df3['helix_len'] == 10]
#     df3['sequences_Helix'].to_csv('data/helix_10.csv', header=True, sep='\t', index=False)
#
#     df3['translated_into_groups'] = df2['sequences_Helix'].apply(lambda x: aa_into_group(x))
#     df3['translated_into_BY'] = df2['sequences_Helix'].apply(lambda x: aa_into_BY(x))
#     df3.to_csv('data/helix_10_groups.csv', header=True, index=False)
#     print('done!')
