import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from ComputationalBiology.data_analysis.all_features_calculator import create_species_df
from ComputationalBiology.data_analysis.kmers_analysis import analyze_kmers
from ComputationalBiology.file_utils import genbank_parser

species_name = 'BS3610'
#species_name = 'UA159'
#species_name = 'AL590842.1_EColi'

overrideFeaturesFile = False
FEATURES_DF_FILE = '../../data/data_outputs/features_' + species_name + '.pickle'

overrideSpeciesParserFile = False
SPECIES_PARSER_FILE = '../../data/data_outputs/species_' + species_name + '.pickle'

overrideKmersDictFile = False
KMERS_DICT_FILE = '../../data/data_outputs/kmers_dict_' + species_name + '.pickle'

if __name__ == '__main__':

    #PARSE
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

    #COMPUTE FEATURES
    if not os.path.exists(FEATURES_DF_FILE) or overrideFeaturesFile:
        species_df = create_species_df(spp)
        species_df.to_pickle(FEATURES_DF_FILE)
    else:
        species_df = pd.read_pickle(FEATURES_DF_FILE)

    #COMPUTE KMERS
    if not os.path.exists(KMERS_DICT_FILE) or overrideKmersDictFile:
        kmers_dict = analyze_kmers(species_df)
        with open(KMERS_DICT_FILE, 'wb') as pickle_file:
            pickle.dump(kmers_dict, file=pickle_file)
    else:
        with open(KMERS_DICT_FILE, 'rb') as pickle_file:
            kmers_dict = pickle.load(file=pickle_file)

    #VISUALIZE (JUNK!)
    kmers_mean = []
    kmers_var = []
    kmers_count = []
    for k in kmers_dict:
        kmers_mean.append(np.mean(kmers_dict[k]))
        kmers_var.append(np.var(kmers_dict[k]))
        kmers_count.append(len(kmers_dict[k]))

    plt.figure(1)
    plt.plot(kmers_mean, 'o')
    plt.show()

    plt.figure(2)
    plt.plot(kmers_var, 'o')
    plt.show()

    plt.figure(3)
    plt.plot(kmers_count, 'o')
    plt.show()

    print('done')

# TODO: operons?
# TODO: assert that locus_tag is unique
#
# # for cds: use codon_start, translation_table
