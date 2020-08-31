import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import NaN
from pandas import DataFrame

from ComputationalBiology.data_analysis.all_features_calculator import GeneFeatures, ProteinFeatures

FEATURES_DF_FILE = '../../data/dataframes/features_BS3610.pickle'


def get_df_by_product(df: pd.DataFrame, product_type: str):
    return df[df['PRODUCT_TYPE'] == product_type]


def plot_hist_feature(df, feature_name: str, out_file=''):

    # filtering out non valid values
    col = df[feature_name]
    x = col[~col.isnull()].values
    ig = sns.distplot(x)
    plt.xlabel(feature_name)
    plt.ylabel('counts')
    plt.title(feature_name)
    plt.show()

    if out_file != '':
        fig.savefig(out_file)


def plot_all_hist_features(df: DataFrame):

    for gene_feature in list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys()):
        assert(gene_feature in df.columns)  # Assuming that the column exists
        print(gene_feature)
        plot_hist_feature(df, gene_feature, out_file='')  # 'path//{}.png'.format(gene_feature.name)


if __name__ == '__main__':
    # Load data:
    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
    df = pd.read_pickle(FEATURES_DF_FILE)

    df[df == None] = NaN

    # Print statistics:
    print(df.describe())

    # print report columns:
    print(df.columns)

    # Get only CDS samples:
    df_cds = get_df_by_product(df, 'CDS')

    # probably a missing data:
    # df_cds[df_cds['HYDROPHOBIC_AA'].isnull()].iloc[0]

    # plot_all_hist_features(df_cds)

    # get only numerical
    table = df_cds.corr()
    print(table)
    print('done')
    exit(1)

    # compute correlation table


    # Plot histogram of lengths:
    # plot_hist_feature(df, feature_name='LENGTH')

    # scatter plot:
    # fig2 = sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)

    # report genes with melting point > 10,000
    # print(df_cds[df_cds['MELTING_POINT'] > 10000][['GENE_NAME', 'PRODUCT_DESCRIPTION']])
    print('Done main_analysis')




