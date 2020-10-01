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
from ComputationalBiology.genetics_project.main_parser_features_calc import species_name, FEATURES_DF_FILE


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
    df_all = pd.read_pickle(FEATURES_DF_FILE)
    df_cds = get_df_by_product(df_all, 'CDS')

    all_kmers = get_all_kmers(6, ValidAlphabet.NT)

    df_kmers_all = df_all[all_kmers]
    df_kmers_cds = df_cds[all_kmers]

    sum_kmers_all = df_kmers_all.sum()
    sum_kmers_cds = df_kmers_cds.sum()

    #Analysis

    #print GC content mean value for species
    gc_mean = df_cds['GC_CONTENT'].mean()
    print(gc_mean)

    # Barplot of sum counts of all kmers - one figure per species
    #fig_bar = plt.bar(list(range(len(all_kmers))), list(sum_kmers_cds[all_kmers].values))
    #plt.show()

    #fig_plot = plt.plot(sum_kmers_cds[all_kmers])
    #plt.show()

    #find frequent kmers
    #frequent_kmers_all = df_kmers_cds[sum_kmers_all[sum_kmers_all > 5000].keys()]
    #frequent_kmers_cds = df_kmers_cds[sum_kmers_cds[sum_kmers_cds > 5000].keys()]

    # print('ALL GENES')
    # print(sum_kmers_all[sum_kmers_all > 5000].sort_values(ascending=False))
    print('CDS GENES')
    frequent_cds_kmers = sum_kmers_cds[sum_kmers_cds > 5000].sort_values(ascending=False)
    print(frequent_cds_kmers)

    # Boxplot of the kmers' distribution over all genes
    #fig_box = sns.boxplot(data=df_cds[frequent_cds_kmers.keys()], showfliers=False)
    # fig_box = plt.boxplot(df_cds[frequent_cds_kmers.keys()].values)
    print(sum_kmers_cds[['AAAAAA', 'CCCCCC', 'GAGAGG', 'TGTGTG']])
    fig_box1 = plt.boxplot(df_cds[['AAAAAA']].values)
    plt.show()

    print('done')
    #frequent_kmers_all = df_kmers_all[sum_kmers_all[sum_kmers_all > 5000].keys()]
    #print(frequent_kmers_all.sum())

    #create features histograms and heatmap
    # plot_all_features_histograms(df_cds)
    # plot_all_features_heatmap(df_cds)


    # scatter plot:
    # fig2 = sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)

    # report genes with melting point > 10,000
    # print(df_cds[df_cds['MELTING_POINT'] > 10000][['GENE_NAME', 'PRODUCT_DESCRIPTION']])
    print('Done main_analysis')




