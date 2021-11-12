# Goal: print to a file all possible translation from the original sequence:
# a. 5 groups
# b. 4 groups
# c. all chemical features groups (i.e. by beans)

# print 4 files: per 2 species (s mutans, s sobrinus) * two types of sequences (transmembrane, all protein)
from tqdm import tqdm

from ComputationalBiology.bio_general.bio_macros import chem_df
from ComputationalBiology.bio_general.bio_utils import get_all_kmers
import pandas as pd

from ComputationalBiology.bio_general.feature_utils import translate_chemical_feature
from ComputationalBiology.biore.biore_utils import aa_into_group, aa_into_hydrophobicity_group


def get_all_kmers_seq(seq, k):
    res = []
    for i in range(len(seq) - k + 1):
        sub_str = seq[i:i + k]
        res.append(sub_str)
    return res


def get_all_kmers_file(file, k):
    # given a file of sequences returns a list of all appearing kmers

    df = pd.read_csv(file, header=0, sep='\t')
    df = df.drop(columns=['five_groups', 'hydrophobic_groups'])
    df['kmers'] = df['seq'].apply(lambda seq: get_all_kmers_seq(seq, k))

    # get the values of kmers which are many lists and merge to one list
    all_kmers = df['kmers'].explode().values

    return all_kmers


if __name__ == '__main__':
    k = 7

    organisms = ['strmu']#, 'strsobrinus']
    for region_type in tqdm(['TRANSMEM']):#, 'ALL_PROT']):
        for organism in tqdm(organisms):

            # A. get all kmers from all sequences
            lst_all_kmers = get_all_kmers_file('data\{}\{}.csv'.format(organism, region_type), k)

            # create a new dataframe with all kemrs (including duplicates) as values
            working_df = pd.DataFrame(data=lst_all_kmers)

            # for renaming column 0:
            working_df['original_seq'] = working_df[0]
            working_df = working_df.drop(columns=[0])

            # start applying computation of all types of translations:
            working_df['five_groups'] = working_df['original_seq'].apply(lambda x: aa_into_group(x))
            working_df['hydrophob_groups'] = working_df['original_seq'].apply(lambda x: aa_into_hydrophobicity_group(x))

            # generate mapping of 5 groups for chemical features
            for feature_name in chem_df.columns:
                working_df[feature_name] = working_df['original_seq'].apply(lambda x: translate_chemical_feature(x, feature_name))
            fname='data/ISF/all_translations/{}_{}_translations'.format(organism, region_type)

            working_df.to_csv(fname + '.txt', '\t', index=False)
            working_df.to_pickle(fname + '.pickle')
            print('done')
    print('done')