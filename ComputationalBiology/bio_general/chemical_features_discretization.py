from ComputationalBiology.bio_general.bio_macros import chem_df
import pandas as pd


def convert_chemical_features_to_discrete(n_bins: int):
    """
    Given the chemical features df, compute for each feature the discrete
    dictionary that describes it (of n_bins bins)
    :return:
    """
    # a. convert all string to numerical values
    chem_df_discrete = chem_df.copy()
    for feature_name in chem_df_discrete.columns:
        chem_df_discrete[feature_name] = chem_df_discrete[feature_name].apply(lambda x: float(x))

    # b. for each feature compute the discrete values, add as a new columns in order to have the same index
    for feature_name in chem_df_discrete.columns:
        chem_df_discrete[feature_name+'_bins'] = pd.cut(chem_df_discrete[feature_name],
                                                        bins=n_bins,
                                                        labels=['A', 'B', 'C', 'D', 'E'],
                                                        right=False)

    return chem_df_discrete


chem_df_discrete = convert_chemical_features_to_discrete(5)

def translate_chemical_feature(seq, feature_name):
    """
    given a sequence of amino acids, translate each letter to the matching bin
    of the requested feature_name
    """
    res = ''
    for c in seq:
        res += str(chem_df_discrete[feature_name+'_bins'].loc[c])

    return res

