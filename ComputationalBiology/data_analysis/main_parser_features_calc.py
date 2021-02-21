import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sklearn
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

from ComputationalBiology.data_analysis.all_features_calculator import create_species_df
from ComputationalBiology.data_analysis.kmers_analysis import create_kmers_df, create_prefix_suffix_agg_df
from ComputationalBiology.file_utils import genbank_parser

species_name = 'BS3610'
#species_name = 'UA159'
#species_name = 'AL590842.1_EColi'
#species_name = 'Bacillus-clausii'

overrideFeaturesFile = True
FEATURES_DF_FILE = '../../data/data_outputs/features_' + species_name + '.pickle'

overrideSpeciesParserFile = True
SPECIES_PARSER_FILE = '../../data/data_outputs/species_' + species_name + '.pickle'

#overrideKmersDictFile = False
#KMERS_DF_FILE = '../../data/data_outputs/kmers_dict_' + species_name + '.pickle'

overridePrefixSuffixDictFile = True
PREFIX_SUFFIX_DF_FILE = '../../data/data_outputs/prefix_suffix_dict_' + species_name + '.pickle'

if __name__ == '__main__':
    # PARSE
    if not os.path.exists(SPECIES_PARSER_FILE) or overrideSpeciesParserFile:
        genbank_file = '../../data/data_inputs/' + species_name + '.gb'
        assert (os.path.exists(genbank_file))  # making sure that the path is valid
        record_gb = genbank_parser.read_genbank_file(genbank_file)
        spp = genbank_parser.init_species(record_gb)
        with open(SPECIES_PARSER_FILE, 'wb') as pickle_file:
            pickle.dump(spp, file=pickle_file)
    else:
        with open(SPECIES_PARSER_FILE, 'rb') as pickle_file:
            spp = pickle.load(file=pickle_file)

    # COMPUTE FEATURES PER GENE
    if not os.path.exists(FEATURES_DF_FILE) or overrideFeaturesFile:
        species_df = create_species_df(spp)
        species_df.to_pickle(FEATURES_DF_FILE)
    else:
        species_df = pd.read_pickle(FEATURES_DF_FILE)

    # COMPUTE NEXT_NT_DF
    if not os.path.exists(PREFIX_SUFFIX_DF_FILE) or overridePrefixSuffixDictFile:
        next_prefix_suffix_agg_dict, prefix_count_agg_dict = create_prefix_suffix_agg_df(species_df, 'CDS', min_gene_length=500, max_gene_length=1000)
        prefix_suffix_tp = (next_prefix_suffix_agg_dict, prefix_count_agg_dict)
        with open(PREFIX_SUFFIX_DF_FILE, 'wb') as pickle_file:
            pickle.dump(prefix_suffix_tp, file=pickle_file)
    else:
        with open(PREFIX_SUFFIX_DF_FILE, 'rb') as pickle_file:
            prefix_suffix_tp = pickle.load(file=pickle_file)

    # # FURTHER ANALYZE STATISTICALLY SIGNIFICANT SEQUENCES FOUND IN PREVIOUS STEP
    # if not os.path.exists(KMERS_DF_FILE) or overrideKmersDictFile:
    #     kmers_df = create_kmers_df(species_df, 'CDS')
    #     with open(KMERS_DF_FILE, 'wb') as pickle_file:
    #         pickle.dump(kmers_df, file=pickle_file)
    # else:
    #     with open(KMERS_DF_FILE, 'rb') as pickle_file:
    #         kmers_df = pickle.load(file=pickle_file)

exit(0)
