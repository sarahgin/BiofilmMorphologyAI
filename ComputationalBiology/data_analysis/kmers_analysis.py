import pandas as pd
import numpy as np

from ComputationalBiology.bio_general.bio_macros import ValidAlphabet
from ComputationalBiology.bio_general.bio_utils import kmers_generator


from ComputationalBiology.data_analysis.all_features_calculator import GeneralFeatures, GeneFeatures


def create_kmers_df(species_df, product_type: str):
    #Step 1: create a dict with hexamers as keys and values
    # are lists of all normalized positions
    kmers_dict = {}
    count = 1
    species_df = species_df[species_df['PRODUCT_TYPE'] == product_type]
    for index, row in species_df.iterrows():
        print(count)
        count += 1
        hex_dict = row[GeneralFeatures.HEXAMER_DICT.name]
        for k in hex_dict:
            hex_norm = np.divide(hex_dict[k], row[GeneFeatures.DNA_LENGTH.name]-len(k)+1)
            if k in kmers_dict:
                kmers_dict[k] = np.append(kmers_dict[k], hex_norm)
            else:
                kmers_dict[k] = hex_norm

    #Step 2: create a list of dictionaries - each dictionary will
    #include {'kmer':'AAAAAA', 'position: [0.2,0.3,0.4]'
    hex_dict_list = []
    for k in kmers_dict:
        hex_dict_list.append({'KMER': k, 'RELATIVE_POSITIONS': kmers_dict[k]})

    return pd.DataFrame(hex_dict_list)
