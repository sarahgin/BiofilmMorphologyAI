import os
import pickle

from ComputationalBiology.data_analysis.all_features_calculator import create_species_df
from ComputationalBiology.file_utils import genbank_parser

species_name = 'BS3610'
# species_name = 'UA159'
# species_name = 'AL590842.1_EColi'
FEATURES_DF_FILE = '../../data/dataframes/features_' + species_name + '.pickle'

if __name__ == '__main__':
    genbank_file = '../../data/GeneticData/' + species_name + '.gb'

    assert(os.path.exists(genbank_file))  # making sure that the path is valid
    record_gb = genbank_parser.read_genbank_file(genbank_file)
    spp = genbank_parser.init_species(record_gb)
    print('done sepcies')
    # TODO: serialize the species
    pickle_file_path = '../../data/SpeciesPickle/species_' + species_name + '.pickle'
    print(pickle_file_path)
    with open(pickle_file_path, 'wb') as pickle_file:
        pickle.dump(spp, file=pickle_file)

    with open(pickle_file_path, 'rb') as pickle_file:
        spp2 = pickle.load(file=pickle_file)

    print('Creating  df of features...')
    df = create_species_df(spp)
    df2 = create_species_df(spp2)  # debug

    import numpy as np
    print(np.all(np.all(df == df2)))
    # Serialize the df
    if not os.path.exists(FEATURES_DF_FILE):
        print('Saving pickle file: {}'.format(FEATURES_DF_FILE))
    else:
        print('Overriding pickle file: {}'.format(FEATURES_DF_FILE))

    df.to_pickle(FEATURES_DF_FILE)

    print('done')


# TODO: operons?
# TODO: assert that locus_tag is unique
#
# # for cds: use codon_start, translation_table

