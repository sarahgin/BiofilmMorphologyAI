import pandas as pd
import numpy as np

from ComputationalBiology.bio_general.bio_macros import ValidAlphabet
from ComputationalBiology.bio_general.bio_utils import kmers_generator


from ComputationalBiology.data_analysis.all_features_calculator import GeneralFeatures, GeneFeatures


def analyze_kmers(species_df):
    #create a dict with hexamers as keys and values
    # are lists of all normalized positions
    kmers_dict = {}
    count = 1
    for index, row in species_df.iterrows():
        print(count)
        count += 1
        hex_dict = row[GeneralFeatures.HEXAMER_DICT.name]
        for k in hex_dict:
            hex_norm = np.divide(hex_dict[k], row[GeneFeatures.DNA_LENGTH.name])
            if k in kmers_dict:
                kmers_dict[k] = np.append(kmers_dict[k], hex_norm)
            else:
                kmers_dict[k] = hex_norm
    return kmers_dict
