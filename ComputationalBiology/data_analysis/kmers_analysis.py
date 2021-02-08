import pandas as pd
import numpy as np

from ComputationalBiology.bio_general.bio_macros import ValidAlphabet
from ComputationalBiology.bio_general.bio_utils import kmers_generator


from ComputationalBiology.data_analysis.all_features_calculator import GeneralFeatures, GeneFeatures


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

        hex_dict = row[GeneralFeatures.HEXAMER_DICT.name]
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

def create_next_nucleotide_df(species_df, product_type: str, min_gene_length=0, max_gene_length=np.inf):
    #Step 1: create a dict with hexamers as keys and values
    # are lists of all normalized positions
    next_nucleotide_dict = {}
    count = 1

    #filter out step
    species_df = species_df[species_df['PRODUCT_TYPE'] == product_type]
    species_df = species_df[species_df['DNA_LENGTH'] >= min_gene_length]
    species_df = species_df[species_df['DNA_LENGTH'] <= max_gene_length]

    for index, row in species_df.iterrows(): #for each Gene
        print(count)
        count += 1

        gene_next_nucleotide_dict = row[GeneralFeatures.HEXAMER_NEXT_NUCLEOTIDE.name]
        for k in gene_next_nucleotide_dict:
            hex_next_nucleotide_dict = gene_next_nucleotide_dict[k]
            if k in next_nucleotide_dict:
                next_nucleotide_dict[k] = next_nucleotide_dict[k] + np.array(list(hex_next_nucleotide_dict.values()))
            else:
                next_nucleotide_dict[k] = np.array(list(hex_next_nucleotide_dict.values()))

    next_nucleotide_df = pd.DataFrame.from_dict(next_nucleotide_dict, orient='index', columns=['A', 'C', 'G', 'T'])

    return next_nucleotide_df


