import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
FEATURES_DF_FILE = '../../data/dataframes/features.pickle'


def get_df_by_product(df: pd.DataFrame, product_type: str):
    return df[df['PRODUCT_TYPE'] == product_type]


if __name__ == '__main__':
    # Load data:
    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
    df = pd.read_pickle(FEATURES_DF_FILE)

    # Print statistics:
    print(df.describe())

    # print report columns:
    print(df.columns)

    # Get only CDS samples:
    df_cds = get_df_by_product(df, 'CDS')
    print('---------')
    print(df_cds.describe())

    # Plot histogram of lengths:
    x = df_cds['LENGTH'].values
    fig = sns.distplot(x)
    plt.xlabel('LENGTH')
    plt.ylabel('counts')
    plt.title('Length histogram of CDS genes')
    plt.show(fig)

    # GC:
    x = df_cds['GC_CONTENT'].values
    fig = sns.distplot(x)
    plt.xlabel('GC_CONTENT')
    plt.ylabel('counts')
    plt.title('GC_CONTENT histogram of CDS genes')
    plt.show(fig)

    # scatter plot:
    # NOTE: I EXPECT A POSITIVE CORRELATION!
    fig2 = sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)
    plt.show(fig2)

    # report genes with melting point > 10,000
    print(df_cds[df_cds['MELTING_POINT'] > 10000][['GENE_NAME', 'PRODUCT_DESCRIPTION']])
    print('Done main_analysis')




