import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from numpy import NaN
from pandas import DataFrame

from ComputationalBiology.bio_general.bio_utils import get_all_kmers
from ComputationalBiology.data_analysis.all_features_calculator import GeneFeatures, ProteinFeatures, kmers_generator, \
    ValidAlphabet
from ComputationalBiology.file_utils.io_utils import create_dir_if_not_exists
from ComputationalBiology.genetics_project.main_parser import species_name, FEATURES_DF_FILE


def get_df_by_product(df: pd.DataFrame, product_type: str):
    return df[df['PRODUCT_TYPE'] == product_type]


def plot_hist_feature(df, feature_name: str, out_file='', show=False):

    # filtering out non valid values
    col = df[feature_name]
    x = col[~col.isnull()].values

    fig = plt.figure()
    col[~col.isnull()].plot.hist(bins=100)
    #fig = sns.distplot(x)
    plt.xlabel(feature_name)
    plt.ylabel('counts')
    plt.title(feature_name)
    if show:
        plt.show()

    if out_file != '':
        fig.savefig(out_file)


def plot_all_features_histograms(df: DataFrame, show=False):

    for gene_feature in list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys()):
        assert(gene_feature in df.columns)  # Assuming that the column exists
        print(gene_feature)
        out_file = '../../results/features_histograms/{}/{}.png'.format(species_name, gene_feature)
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
    out_file = '../../results/features_correlations/{}/{}.png'.format(species_name, 'corr')
    create_dir_if_not_exists(out_file)
    fig.savefig(out_file)


if __name__ == '__main__':
    # Load data and prepare df_cds:
    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
    df = pd.read_pickle(FEATURES_DF_FILE)
    print(df.describe())
    print(df.columns)
    df_cds = get_df_by_product(df, 'CDS')
    all_kmers = get_all_kmers(6, ValidAlphabet.NT)

    #print GC content mean value for species
    gc_mean = df_cds['GC_CONTENT'].mean()
    print(gc_mean)

    # Barplot of sum counts of all kmers - one figure per species
    # fig = plt.bar(list(range(len(all_kmers))), s[all_kmers])
    # plt.show()

    # Boxplot of the kmers' distribution over all genes
    #sns.boxplot(data=df_cds[all_kmers])

    #find frequent kmers
    s = df_cds.sum() # create sum per column
    df_kmers = df_cds[all_kmers]  # get only kmers columns
    sum_kmers = df_kmers.sum() # vector of 4096 elements
    frequent_kmers = df_kmers[sum_kmers[sum_kmers > 3000].keys()]
    print(frequent_kmers.columns)

    #create features histograms and heatmap
    # plot_all_features_histograms(df_cds)
    # plot_all_features_heatmap(df_cds)


    # scatter plot:
    # fig2 = sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)

    # report genes with melting point > 10,000
    # print(df_cds[df_cds['MELTING_POINT'] > 10000][['GENE_NAME', 'PRODUCT_DESCRIPTION']])
    print('Done main_analysis')




