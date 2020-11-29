import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sklearn
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

from ComputationalBiology.data_analysis.all_features_calculator import create_species_df
from ComputationalBiology.data_analysis.kmers_analysis import create_kmers_df
from ComputationalBiology.file_utils import genbank_parser

species_name = 'BS3610'
# species_name = 'UA159'
# species_name = 'AL590842.1_EColi'

overrideFeaturesFile = False
FEATURES_DF_FILE = '../../data/data_outputs/features_' + species_name + '.pickle'

overrideSpeciesParserFile = False
SPECIES_PARSER_FILE = '../../data/data_outputs/species_' + species_name + '.pickle'

overrideKmersDictFile = False
KMERS_DF_FILE = '../../data/data_outputs/kmers_dict_' + species_name + '.pickle'

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

    # COMPUTE FEATURES
    if not os.path.exists(FEATURES_DF_FILE) or overrideFeaturesFile:
        species_df = create_species_df(spp)
        species_df.to_pickle(FEATURES_DF_FILE)
    else:
        species_df = pd.read_pickle(FEATURES_DF_FILE)

    # COMPUTE KMERS
    if not os.path.exists(KMERS_DF_FILE) or overrideKmersDictFile:
        kmers_df = create_kmers_df(species_df, 'CDS', 500, 2000)
        with open(KMERS_DF_FILE, 'wb') as pickle_file:
            pickle.dump(kmers_df, file=pickle_file)
    else:
        with open(KMERS_DF_FILE, 'rb') as pickle_file:
            kmers_df = pickle.load(file=pickle_file)

    #KMeans analysis
    kmers_df['POSITIONAL_BINS_DISTRIBUTION'] \
        = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: np.histogram(x, bins=10)[0] / np.sum(x))

    array_to_cluster_list = kmers_df['POSITIONAL_BINS_DISTRIBUTION'].to_list()
    #kmeans = KMeans(n_clusters=2, random_state=0).fit(array_to_cluster.to_list())

    #Elbow method
    distortions = []
    silhouette_scores = []
    K = range(2, 20)
    for k in K:
        print(k)
        kmeanModel = KMeans(n_clusters=k, random_state=0).fit(array_to_cluster_list)

        score = sklearn.metrics.silhouette_score(array_to_cluster_list, kmeanModel.labels_)

        silhouette_scores.append(score)
        distortions.append(sum(np.min(cdist(array_to_cluster_list, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) /
                           len(array_to_cluster_list))

    # Plot the elbow
    plt.plot(K, distortions, 'bx-')
    plt.plot(K, silhouette_scores, 'rx-')
    plt.xlabel('k')
    plt.ylabel('Distortion')
    plt.title('The Elbow Method showing the optimal k')
    plt.show()


    print('done')

# TODO: operons?
# TODO: assert that locus_tag is unique
#
# # for cds: use codon_start, translation_table
