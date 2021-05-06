import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
import re

from ComputationalBiology.data_analysis.all_features_calculator import GeneFeatures, ProteinFeatures
from ComputationalBiology.file_utils.io_utils import create_dir_if_not_exists
from ComputationalBiology.data_analysis.main_parser_features_calc import species_name, FEATURES_DF_FILE

from sklearn.decomposition import PCA
# compute PCA on df_joined
from sklearn.preprocessing import StandardScaler


def plot_hist_feature(df, feature_name: str, out_file='', show=False):
    # filtering out non valid values
    col = df[feature_name]
    x = col[~col.isnull()].values

    fig = plt.figure()
    col[~col.isnull()].plot.hist(bins=100)
    # fig = sns.distplot(x)
    plt.xlabel(feature_name)
    plt.ylabel('counts')
    plt.title(feature_name)
    if show:
        plt.show()

    if out_file != '':
        fig.savefig(out_file)


def plot_all_features_histograms(df: DataFrame, show=False, suffix=''):
    for gene_feature in list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys()):
        print(gene_feature)
        assert (gene_feature in df.columns)  # Assuming that the column exists
        out_file = '../../data/data_graphs/features_histograms/{}/{}.png'.format(species_name + suffix, gene_feature)
        create_dir_if_not_exists(out_file)
        plot_hist_feature(df, gene_feature, out_file=out_file, show=show)


def plot_all_features_heatmap(df, show=False, suffix=''):
    # plot heatmap of correlations
    df = df.replace(to_replace='None', value=np.nan).dropna()
    table = df.corr()
    fig = plt.figure(figsize=(8.0, 5.0))
    sns.heatmap(table, annot=True)
    if show:
        plt.show()
    out_file = '../../data/data_graphs/features_correlations/{}/{}.png'.format(species_name + suffix, 'corr')
    create_dir_if_not_exists(out_file)
    fig.savefig(out_file)


def plot_pca(df_joined, target_column):#='is_gene_of_interest'):

    features =[ \
        # 'GC_CONTENT',
     # 'DNA_LENGTH',
     # 'HYDROPHOBIC_AA',
     # 'HYDROPHILIC_AA',
     # 'POLAR_AA',
     'AROMATIC_AA',
     'POSITIVE_AA',
     # 'NEGATIVE_AA',
     'NONPOLAR_AA',
     # 'AA_LENGTH',
     'H1',
     # 'H2',
     # 'H3',
     # 'V',
     # 'P1',
     'P2',
     'SASA',
     # 'NCI',
     'MASS'
     ]

    # Separating out the features
    df_joined = df_joined.dropna(0) # NOTE: in BS3610 we had 58 genes of NAN in some of the features
    x = df_joined.loc[:, features].values

    # Separating out the target
    y = df_joined.loc[:, [target_column]].values

    # Standardizing the features
    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x, )
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])

    finalDf = pd.concat([principalDf, df_joined[[target_column]]], axis=1)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    # targets = [False, True]
    # targets = list(set(finalDf[target_column].values))

    # remove the most frequent labels and add to the end
    from collections import Counter
    label_counter_sorted = Counter(finalDf[target_column].values).most_common()
    targets = [x[0] for x in label_counter_sorted]
    print(Counter(finalDf[target_column].values))

    colors = ['thistle', 'olive', 'skyblue', 'firebrick', 'gold', 'darkred']
    for target, color, in zip(targets, colors):
        indicesToKeep = finalDf[target_column] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'],
                   finalDf.loc[indicesToKeep, 'principal component 2'],
                   s=50, c=color)
    ax.legend(targets)
    ax.grid()
    plt.show()
    print('done')


#FOR POSTER GRAPHS ONLY
if __name__ == '__main__':
    # FEATURES_DF_FILE = '../../data/data_outputs/features_BS168.pickle'
    # FEATURES_DF_FILE
    df_all_BS168 = pd.read_pickle(FEATURES_DF_FILE)
    print(FEATURES_DF_FILE)
    df_cds_BS168 = df_all_BS168[df_all_BS168['PRODUCT_TYPE'] == 'CDS'].copy()
    # s1 = set(df_cds_BS168['GENE_NAME'])
    # change to lower in order to join with the uniprot table
    df_cds_BS168['GENE_NAME'] = df_cds_BS168['GENE_NAME'].apply(lambda x: x.lower())

    use_uniprot = False
    if use_uniprot:
        # UNIPROT DF:
        UNIPROT_FILE = '../../data/data_inputs/bs168_uniprot.tab'
        df_UP = pd.read_csv(UNIPROT_FILE, sep='\t', header=0)
        df_UP = df_UP.fillna('')
        df_UP['isSecreted'] = df_UP['Subcellular location [CC]'].apply(lambda x: re.search('Secreted', x) is not None)
        df_UP['isExtracellular'] = df_UP['Topological domain'].apply(lambda x: re.search('Extracellular', x) is not None)
        # Note: we fetched only the 1st name in the list
        df_UP['GENE_NAME'] = df_UP['Gene names'].apply(lambda x: x.lower().split()[0] if x != '' else 'NA_GENE_NAME_UNIPROT')

        df_joined = pd.merge(df_cds_BS168, df_UP, on="GENE_NAME")

        df_joined['is_gene_of_interest'] = df_joined['isSecreted'] | df_joined['isExtracellular']
        df_genes_interest = df_joined[df_joined['isSecreted'] | df_joined['isExtracellular']]
        df_other_genes = df_joined[~df_joined['isSecreted'] & ~df_joined['isExtracellular']]
    else:
        df_joined = df_cds_BS168
        df_joined = df_joined.reset_index()

    # plot_all_features_histograms(df_genes_interest, suffix='_Secreted_Extracellular')
    # plot_all_features_histograms(df_other_genes, suffix='_Other')
    # plot_all_features_heatmap(df_genes_interest, suffix='_Secreted_Extracellular')
    # plot_all_features_heatmap(df_other_genes, suffix='_Other')

    # set all labels to true just in order to have a single color

    df_joined['is_gene_of_interest'] = True

    # adding new labels from:
    # http://subtiwiki.uni-goettingen.de/wiki/index.php/Biofilm_formation#Key_genes_and_operons_involved_in_biofilm_formation
    df_joined['is_gene_of_biofilm'] = 'other'

    dict_biofilm = {
                    'polysaccharide': ['epsA', 'epsB', 'epsC', 'epsD', 'epsE', 'epsF', 'epsG',
                                       'epsH', 'epsI', 'epsJ', 'epsK', 'epsL', 'epsM', 'epsN', 'epsO', 'galE'],
                    'amyloid': ['tapA', 'sipW', 'tasA'],
                    'regulation': ['SlrR', 'SlrA', 'SinR', 'SinI', 'KinC', 'KinD', 'Spo0A', 'PtkA', 'TkmA', 'PtpZ',
                     'DegU', 'DegQ', 'YmdB', 'FtsH', 'Veg', 'MstX', 'YugO'],
                    'formation': ['AmpS', 'FloT', 'LuxS',  'RemA', 'RemB', 'Rny', 'Sfp/1', 'Sfp/2',
                                  'SpeA', 'SpeD', 'SwrAA', 'YisP', 'YlbF', 'YmcA', 'YvcA', 'YwcC', 'YxaB'],
                    'pellicle': ['Hag', 'FlgE', 'FliF', 'MotA', 'SigD', 'CheA',
                     'CheY', 'CheD', 'CheV', 'HemAT']

                    }
    # # epsA - epsB - epsC - epsD - epsE - epsF - epsG - epsH - epsI - epsJ - epsK - epsL - epsM - epsN - epsO
    # genes_biofilm = ['galE', 'tapA', 'sipW', 'tasA', 'bslA']
    dict_biofilm_lower = {}
    for k in dict_biofilm.keys():
        dict_biofilm_lower[k] = [x.lower() for x in dict_biofilm[k]]
    # genes_biofilm = [x.lower() for x in genes_biofilm]
    # print('len(genes_biofilm):', len(genes_biofilm))
    count_found = 0
    for k_group in dict_biofilm_lower.keys():
        for gene in dict_biofilm_lower[k_group]:
            count_found += sum(df_joined['GENE_NAME'] == gene)
            df_joined.loc[df_joined['GENE_NAME'] == gene, ['is_gene_of_biofilm']] = k_group
    print(count_found)
    # print(sum(df_joined['is_gene_of_biofilm']))
    plot_pca(df_joined, target_column='is_gene_of_biofilm') #'is_gene_of_interest')

    print('done')
    # sres = set(s1_lower).intersection(set(s2_lower))
    # print('done')

# ORIGINAL MAIN:
#if __name__ == '__main__':
#    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
#    df_all = pd.read_pickle(FEATURES_DF_FILE)
#    df_cds = df_all[df_all['PRODUCT_TYPE'] == 'CDS']

#    # create features histograms and heatmap
#    plot_all_features_histograms(df_cds)
#    plot_all_features_heatmap(df_cds)

#    # scatter plot:
#    #sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)

#    print('Done main_analysis')
