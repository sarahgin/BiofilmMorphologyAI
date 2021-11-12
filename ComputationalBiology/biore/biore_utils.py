from ComputationalBiology.biore.biore_macros import AA_GROUP, group_to_aa_dict, aa_group_list, aa_to_group_dict, \
    aa_to_BY_dict, AA, NT, ValidAlphabet, aa_to_hydrophobicity_group_dict
from scipy.spatial import distance
import pandas as pd
import numpy as np
from scipy.stats import entropy


def translate(original_str: str):
    res = ''
    i = 0
    while i < len(original_str):
        c = original_str[i]
        if c == "\\":
            res += original_str[i + 1]
            i += 1
        elif c in aa_group_list:
            res += group_into_aa_or(c)
        else:
            res += c
        i += 1
    return res


def group_into_aa_or(group_letter: str):
    # TODO: handle case of invalid group letter
    relevant_aa = group_to_aa_dict[group_letter]
    res = '('
    for i in range(len(relevant_aa) - 1):
        res += relevant_aa[i] + '|'
    res += relevant_aa[i + 1] + ')'
    return res


# simple translations single-letter aa into amino acid group letter
def aa_into_group(seq):
    res = ''
    for c in seq:
        res += aa_to_group_dict[c]
    return res


def aa_into_hydrophobicity_group(seq):
    res = ''
    for c in seq:
        res += aa_to_hydrophobicity_group_dict[c]
    return res


def aa_into_group_ignore_dots(seq):
    res = ''
    for c in seq:
        if c == '.':
            res += '.'
        else:
            res += group_into_aa_or(aa_to_group_dict[c])
    return res


def aa_into_BY(seq):
    res = ''
    for c in seq:
        res += aa_to_BY_dict[c]
    return res


def group_into_aa(seq):
    res = ''
    for c in seq:
        res += group_to_aa_dict[c]
    return res


def get_hamming(seq1, seq2):
    assert len(seq1) == len(seq2)
    return distance.hamming(list(seq1), list(seq2))


def get_positional_probability_matrix(sequences, seq_len, alphabet):

    alphabet_enum = AA if alphabet == ValidAlphabet.AA else NT

    df = pd.DataFrame(sequences, columns=['Seq'])
    for p in range(seq_len):
        df['position_' + str(p)] = df['Seq'].apply(lambda x: x[p])

    df_probability_matrix = pd.DataFrame(np.zeros((seq_len, len(alphabet_enum))), columns=[e.value for e in alphabet_enum])

    for p in range(seq_len):
        for g in [e.value for e in alphabet_enum]:
            df_probability_matrix.loc[p, g] = (sum(df['position_' + str(p)] == g))

    seq_num = df_probability_matrix.sum(axis=1).iloc[0]
    df_probability_matrix = df_probability_matrix.applymap(lambda i: i / seq_num)
    return df_probability_matrix

def get_entropy_vector(df_probability_matrix):
    df_probability_matrix['entropy_value'] = \
        df_probability_matrix.apply(lambda row: entropy(row, base=2), axis=1)
    return df_probability_matrix['entropy_value'].values

#TODO: fix code duplication
def create_regex_all_or_none(df_probability_matrix):
    regex = ''
    for pos in range(0, len(df_probability_matrix)):
        relevant_aa = df_probability_matrix.columns[df_probability_matrix.iloc[pos] != 0]
        if len(relevant_aa) == 1:
            regex += relevant_aa[0]
        else: #more than one amino acid found at position
            regex += '.'
    return regex

def create_regex_or(df_probability_matrix):
    regex = ''
    for pos in range(0, len(df_probability_matrix)):
        relevant_aa = df_probability_matrix.columns[df_probability_matrix.iloc[pos] != 0]
        if len(relevant_aa) == 1:
            regex += relevant_aa[0]
        else: #more than one amino acid found at position
            regex += '('
            regex += '|'.join(relevant_aa)
            regex += ')'

    return regex

if __name__ == '__main__':
    #sequences = ['ATGC', 'AATT', 'ATGG', 'ATGA']

    sequences = ['QLYFFSSLVCTECVLLASMAY',
    'QLYFFHFLGSTECFLYTVMSY',
    'QMYFFFFFGVAECFLLATMAY',
    'QMFFFLFFGATECCLLAAMAY',
    'QMFFFVFLGAAECFLLSSMAY',
    'QMFFFVTLGSTDCFLLAIMAY',
    'QLYFFHFLGSTQCFLYTLMAY',
    'QLYFFLGLGCTECVLLAVMAY',
    'QMYFFFFFGAAECCLLATMAY',
    'QIYFFHSLGATECYLLTAMAY',
    'QLYFFIALACTECVLLAVMAY',
    'QLYFFIALMCTECVLLAAMAY',
    'QLYFFVFLGATECFLLAFMAY',
    'QLFIFTFLGATECFLLAAMAY']


    m = get_positional_probability_matrix(sequences, 21, ValidAlphabet.AA)
    print(m)
    reg = create_regex_all_or_none(m)
    print(reg)
    v = get_entropy_vector(m)
    print(v)
