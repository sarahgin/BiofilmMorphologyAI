# goal: get human1.tab -->  fetch all helix
# save in dataframe
import pandas as pd
import re

from ComputationalBiology.biore.biore_utils import aa_into_group, aa_into_BY


def parse_file(file_path):
    """
    The function parses the given tab-separated file and returns a dataframe object
    :return:
    """
    df = pd.read_csv(file_path, sep='\t', header=0)
    print(df.head())
    return df


def get_subsequences(motif_info: str, sequence: str):
    """
    Find all subsequences in sequence according to the ranges defines in motif_info
    :param motif_info: String describing motif(s) locations within sequence
    :param sequence: query sequence
    :return: a list of subsequences
    """
    res = []
    # find all matches of (HELIX (\d+)..(\d+))# TODO: make general fo beta
    all_matches = re.findall('HELIX (\d+)\\.\\.(\d+)', motif_info)
    assert(all_matches is not None)

    # go over the matches
    for start, end in all_matches:
        # print(int(start)- int(end))
        res.append(sequence[int(start)-1:int(end)+1-1])

    return res


def motif_sequences_to_file(df: pd.DataFrame, col_name: str):
    """
    Given a dataframe and a column name (e.g. helix), returns a list
    all of the subtrings that match the motif
    :param col_name:
    :return: a list
    """
      # add a column to check is there are matches of the required column
    # df.fillna(None)
    df['has_' + col_name] = df[col_name].apply(lambda x: type(x) is str)
    df2 = df.copy()

    # remove unneeded rows
    df2 = df2[df2['has_' + col_name]]

    # fetch all subsequences
    df2['sequences_' + col_name] = df2.apply(lambda row: get_subsequences(row[col_name], row['Sequence']), axis=1)

    df3 = df2['sequences_' + col_name].copy()
    # df_all_helices = df3.apply(pd.Series).stack().reset_index(drop=True)
    df_all_helices = df3.explode()

    return df_all_helices


if __name__ == '__main__':
    # parsing file to create all helices:
    if False:
        file_path = 'data/human1.tab'
        df = parse_file(file_path)

        df_all_helices = motif_sequences_to_file(df=df, col_name='Helix')

        df_all_helices.to_csv('data/all_helix.csv', header=True, sep='\t', index=False)

    #read helix data:
    print('reading...')
    df2 = pd.read_csv('data/all_helix.csv')
    df2['helix_len'] = df2['sequences_Helix'].apply(lambda x: len(x))
    df2.describe()

    # get only 10-mers
    df3 = df2.copy()
    df3 = df3[df3['helix_len'] == 10]
    df3['sequences_Helix'].to_csv('data/helix_10.csv', header=True, sep='\t', index=False)

    df3['translated_into_groups'] = df2['sequences_Helix'].apply(lambda x: aa_into_group(x))
    df3['translated_into_BY'] = df2['sequences_Helix'].apply(lambda x: aa_into_BY(x))
    df3.to_csv('data/helix_10_groups.csv', header=True, index=False)
    print('done')