import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame

from ComputationalBiology.bio_general.bio_utils import get_all_kmers
from ComputationalBiology.data_analysis.all_features_calculator import GeneFeatures, ProteinFeatures, ValidAlphabet
from ComputationalBiology.file_utils.io_utils import create_dir_if_not_exists
from ComputationalBiology.data_analysis.main_parser_features_calc import species_name, FEATURES_DF_FILE, KMERS_DF_FILE


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
        assert (gene_feature in df.columns)  # Assuming that the column exists
        print(gene_feature)
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


if __name__ == '__main__':
    # Load data and prepare df_cds:
    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
    df_all = pd.read_pickle(FEATURES_DF_FILE)
    df_cds = df_all[df_all['PRODUCT_TYPE'] == 'CDS']

    # Load kmer df file into df_kmers
    kmers_dict_list = pd.read_pickle(KMERS_DF_FILE)
    kmers_df = pd.DataFrame(kmers_dict_list)
    kmers_df['MEAN'] = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: np.mean(x))
    kmers_df['VAR'] = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: np.var(x))
    kmers_df['COUNT'] = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: len(x))

    # top frequent kmers position distribution
    frequent_kmers = kmers_df[kmers_df['COUNT'] > 5000].sort_values(by='COUNT', ascending=False)
    plt.boxplot(frequent_kmers.iloc[0:100]['RELATIVE_POSITIONS'])
    plt.show()

    # boxplots of the three columns
    kmers_df.boxplot(column=['MEAN'])
    plt.title('MEAN')
    plt.show()
    kmers_df.boxplot(column=['VAR'])
    plt.title('VAR')
    plt.show()
    kmers_df.boxplot(column=['COUNT'])
    plt.title('COUNT')
    plt.show()

    # create features histograms and heatmap
    plot_all_features_histograms(df_cds)
    plot_all_features_heatmap(df_cds)

    # scatter plot:
    # sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)

    print('Done main_analysis')
