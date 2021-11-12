import os
import pickle
import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# import sklearn
# from sklearn.cluster import KMeans
# from scipy.spatial.distance import cdist

# from ComputationalBiology.data_analysis.all_features_calculator import create_species_df
from ComputationalBiology.bio_general.bio_macros import PREFIX_LENGTH_MIN, SUFFIX_LENGTH_MIN, PREFIX_LENGTH_MAX, \
    SUFFIX_LENGTH_MAX
from ComputationalBiology.data_analysis.all_features_calculator import create_species_df
from ComputationalBiology.data_analysis.kmers_analysis import create_kmers_df, create_prefix_suffix_agg_df
from ComputationalBiology.file_utils import genbank_parser

#species_name = 'actinomyces_israelii'
#species_name = 'bacillus_clausii'
species_name = 'bacillus_subtilis'
#species_name = 'escherichia_coli'
##species_name = 'staph_aureus'
#species_name = 'strep_mitis'
# species_name = 'strep_mutans'
##species_name = 'strep_salivarius'
# species_name = 'strep_sanguinis'
# species_name = 'strep_sobrinus'
##species_name = 'helicobacter_pylori'

overrideSpeciesParserFile = True
SPECIES_PARSER_FILE = '../../data/data_outputs/species_' + species_name + '.pickle'

overrideFeaturesFile = True
FEATURES_DF_FILE = '../../data/data_outputs/features_' + species_name + '.pickle'

#overrideKmersDictFile = False
#KMERS_DF_FILE = '../../data/data_outputs/kmers_dict_' + species_name + '.pickle'

#overridePrefixSuffixDictFile = True
#PREFIX_SUFFIX_DF_FILE = '../../data/data_outputs/prefix_suffix_dict_' + species_name + '.pickle'

if __name__ == '__main__':
    # PARSE
    if not os.path.exists(SPECIES_PARSER_FILE) or overrideSpeciesParserFile:
        genbank_file = '../../data/data_inputs/GenBank/' + species_name + '.gb'
        assert (os.path.exists(genbank_file))  # making sure that the path is valid
        record_gb = genbank_parser.read_genbank_file(genbank_file)
        spp = genbank_parser.init_species(record_gb)
        with open(SPECIES_PARSER_FILE, 'wb') as pickle_file:
            pickle.dump(spp, file=pickle_file)
    else:
        with open(SPECIES_PARSER_FILE, 'rb') as pickle_file:
            spp = pickle.load(file=pickle_file)

    # compute features file
    if not os.path.exists(FEATURES_DF_FILE) or overrideFeaturesFile:
        species_df = create_species_df(spp)
        species_df.to_pickle(FEATURES_DF_FILE)
    else:
        species_df = pd.read_pickle(FEATURES_DF_FILE)

print('all done :)')