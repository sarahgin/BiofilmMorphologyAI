import os.path

import pandas as pd
from itertools import product

from BioinformaticsLab.ComputationalBiology.bio_general.bio_macros import ValidAlphabet, alphabets_dict


def kmers_generator(k: int, alphabet: ValidAlphabet):
    """
    This is a generator, generates all DNA k-mers
    :param alphabet:
    :return: the next k-mer (using yield)
    :param k:
    """
    gen = product(alphabets_dict[alphabet], repeat=k)
    for x in gen:
        yield ''.join(x)


def get_all_kmers(k, alphabet):
    all_kmers = []
    for kmer in kmers_generator(k, alphabet):
        all_kmers.append(kmer)
    return all_kmers


def merge_add_dicts(dict1, dict2):
    merged_dict = {key: dict1.get(key, 0) + dict2.get(key, 0)
                   for key in set(dict1) | set(dict2)}
    return merged_dict


def get_correlation_matrix(df):
    # Note: if any of the features is constant value, pearson correlation
    # is not defined, thus its heatmap row will be empty

    columns = ['GC_CONTENT', 'GENE_LENGTH',  # DNA features
               'HAS_TATAAT', 'HAS_TTTATT',  # DNA motif (i.e., bioregex)
               'HAS_PPKL',  # Protein motif (i.e., bioregex)
               # Protein features:
                'PROTEIN_LENGTH', 'HYDROPHOBIC_AA', 'HYDROPHILIC_AA', 'POLAR_AA',
               'AROMATIC_AA', 'POSITIVE_AA', 'NEGATIVE_AA', 'NONPOLAR_AA',
               'H1', 'H2', 'H3', 'V', 'P1', 'P2', 'SASA', 'NCI',
               'MASS', 'PKA_COOH', 'PKA_NH', 'PI']

    # fetch the required columns:
    df = df[columns]
    return df.corr()
