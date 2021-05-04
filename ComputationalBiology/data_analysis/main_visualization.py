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


def plot_all_features_histograms(df: DataFrame, show=False):
    for gene_feature in list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys()):
        print(gene_feature)
        assert (gene_feature in df.columns)  # Assuming that the column exists
        out_file = '../../data/data_graphs/features_histograms/{}/{}.png'.format(species_name, gene_feature)
        create_dir_if_not_exists(out_file)
        plot_hist_feature(df, gene_feature, out_file=out_file, show=show)


def plot_all_features_heatmap(df, show=False):
    # plot heatmap of correlations
    df = df.replace(to_replace='None', value=np.nan).dropna()
    table = df.corr()
    fig = plt.figure(figsize=(8.0, 5.0))
    sns.heatmap(table, annot=True)
    if show:
        plt.show()
    out_file = '../../data/data_graphs/features_correlations/{}/{}.png'.format(species_name, 'corr')
    create_dir_if_not_exists(out_file)
    fig.savefig(out_file)

#FOR POSTER GRAPHS ONLY
if __name__ == '__main__':
    FEATURES_DF_FILE_BS168 = '../../data/data_outputs/features_BS168.pickle'
    df_all_BS168 = pd.read_pickle(FEATURES_DF_FILE_BS168)
    df_cds_BS168 = df_all_BS168[df_all_BS168['PRODUCT_TYPE'] == 'CDS']
    s1 = set(df_cds_BS168['GENE_NAME'])
    s1_lower = []
    for name in s1:
        s1_lower.append(name.lower())

    UNIPROT_FILE = '../../data/data_inputs/bs168_uniprot.tab'
    df_UP = pd.read_csv(UNIPROT_FILE, sep='\t', header=0)
    df_UP = df_UP.fillna('')
    df_UP['isSecreted'] = df_UP['Subcellular location [CC]'].apply(lambda x: re.search('Secreted', x) is not None)
    df_UP['isExtracellular'] = df_UP['Topological domain'].apply(lambda x: re.search('Extracellular', x) is not None)
    df_UP['singleGeneName'] = df_UP['Gene names'].apply(lambda x: x.lower().split()[0])


    # TODO:
    s2 = set(df_UP['Gene names'])
    s2_lower = []
    for name in s2:
        if isinstance(name, str):
            gene_name = name.lower()
            gene_name = gene_name.split()[0]
            s2_lower.append(gene_name)

    sres = set(s1_lower).intersection(set(s2_lower))
    print('done')

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
