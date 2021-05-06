import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
import re

from ComputationalBiology.data_analysis.all_features_calculator import GeneFeatures, ProteinFeatures
from ComputationalBiology.file_utils.io_utils import create_dir_if_not_exists
from ComputationalBiology.data_analysis.main_parser_features_calc import species_name, FEATURES_DF_FILE


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

#FOR POSTER GRAPHS ONLY
if __name__ == '__main__':
    FEATURES_DF_FILE_BS168 = '../../data/data_outputs/features_BS168.pickle'
    df_all_BS168 = pd.read_pickle(FEATURES_DF_FILE_BS168)
    df_cds_BS168 = df_all_BS168[df_all_BS168['PRODUCT_TYPE'] == 'CDS'].copy()
    # s1 = set(df_cds_BS168['GENE_NAME'])
    # change to lower in order to join with the uniprot table
    df_cds_BS168['GENE_NAME'] = df_cds_BS168['GENE_NAME'].apply(lambda x: x.lower())

    # UNIPROT DF:
    UNIPROT_FILE = '../../data/data_inputs/bs168_uniprot.tab'
    df_UP = pd.read_csv(UNIPROT_FILE, sep='\t', header=0)
    df_UP = df_UP.fillna('')
    df_UP['isSecreted'] = df_UP['Subcellular location [CC]'].apply(lambda x: re.search('Secreted', x) is not None)
    df_UP['isExtracellular'] = df_UP['Topological domain'].apply(lambda x: re.search('Extracellular', x) is not None)
    # Note: we fetched only the 1st name in the list
    df_UP['GENE_NAME'] = df_UP['Gene names'].apply(lambda x: x.lower().split()[0] if x != '' else 'NA_GENE_NAME_UNIPROT')

    df_joined = pd.merge(df_cds_BS168, df_UP, on="GENE_NAME")

    df_genes_interest = df_joined[df_joined['isSecreted'] | df_joined['isExtracellular']]
    df_other_genes = df_joined[~df_joined['isSecreted'] & ~df_joined['isExtracellular']]

    plot_all_features_histograms(df_genes_interest, suffix='_Secreted_Extracellular')
    plot_all_features_histograms(df_other_genes, suffix='_Other')
    plot_all_features_heatmap(df_genes_interest, suffix='_Secreted_Extracellular')
    plot_all_features_heatmap(df_other_genes, suffix='_Other')

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
