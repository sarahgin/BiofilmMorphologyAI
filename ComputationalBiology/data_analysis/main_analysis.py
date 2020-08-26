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
    sns.distplot(x)
    plt.show()

    print('Done main_analysis')




