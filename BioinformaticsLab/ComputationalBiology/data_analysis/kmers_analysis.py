import numpy as np
import pandas as pd

from BioinformaticsLab.ComputationalBiology.bio_general.bio_utils import merge_add_dicts
from BioinformaticsLab.ComputationalBiology.data_analysis.all_features_calculator import DNAFeatures


def create_kmers_df(species_df, product_type: str, min_gene_length=0, max_gene_length=np.inf):
    #Step 1: create a dict with hexamers as keys and values
    # are lists of all normalized positions
    kmers_dict = {}
    count = 1

    #filter out step
    species_df = species_df[species_df['PRODUCT_TYPE'] == product_type]
    species_df = species_df[species_df['DNA_LENGTH'] >= min_gene_length]
    species_df = species_df[species_df['DNA_LENGTH'] <= max_gene_length]

    for index, row in species_df.iterrows(): #for each Gene
        print(count)
        count += 1

        hex_dict = row[DNAFeatures.HEXAMER_DICT.name]
        for k in hex_dict:
            hex_norm = hex_dict[k]
            if k in kmers_dict:
                kmers_dict[k] = np.append(kmers_dict[k], hex_norm)
            else:
                kmers_dict[k] = hex_norm

    #Step 2: create a list of dictionaries - each dictionary will
    #include {'kmer':'AAAAAA', 'position: [0.2,0.3,0.4]'
    hex_dict_list = []
    for k in kmers_dict:
        hex_dict_list.append({'KMER': k, 'ABSOLUTE_POSITIONS': kmers_dict[k]})

    return pd.DataFrame(hex_dict_list)


def create_prefix_suffix_agg_df(species_df, product_type='CDS', min_gene_length=0, max_gene_length=np.inf):
    #Step 1: create a dict with hexamers as keys and values
    # are lists of all normalized positions
    next_prefix_suffix_agg_dict = {}
    prefix_count_agg_dict = {}
    count = 0

    #filter out step
    if 'PRODUCT_TYPE' in species_df.columns:
        species_df = species_df[species_df['PRODUCT_TYPE'] == product_type]
    if 'DNA_LENGTH' in species_df.columns:
        species_df = species_df[species_df['DNA_LENGTH'] >= min_gene_length]
        species_df = species_df[species_df['DNA_LENGTH'] <= max_gene_length]

    print('Total number of genes to aggregate: ', len(species_df))
    for index, row in species_df.iterrows(): #for each Gene
        count += 1
        if count % 100 == 0:
            print("create_prefix_suffix_agg_df: gene num-", str(count))

        gene_prefix_suffix_dict, prefix_count = row[KmerFeatures.PREFIX_SUFFIX_DICT.name]
        for k in gene_prefix_suffix_dict:

            if k in prefix_count_agg_dict:
                prefix_count_agg_dict[k] += prefix_count[k]
            else:
                prefix_count_agg_dict[k] = prefix_count[k]

            prefix_suffix_dict = gene_prefix_suffix_dict[k]
            if k in next_prefix_suffix_agg_dict:
                next_prefix_suffix_agg_dict[k] = merge_add_dicts(next_prefix_suffix_agg_dict[k], prefix_suffix_dict)
            else:
                next_prefix_suffix_agg_dict[k] = prefix_suffix_dict

    return next_prefix_suffix_agg_dict, prefix_count_agg_dict


