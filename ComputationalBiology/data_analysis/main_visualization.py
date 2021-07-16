import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
import re

from ComputationalBiology.bio_general.bio_macros import species_names
from ComputationalBiology.data_analysis.all_features_calculator import GeneFeatures, ProteinFeatures
from ComputationalBiology.file_utils.io_utils import create_dir_if_not_exists
from ComputationalBiology.data_analysis.main_parser_features_calc import species_name, FEATURES_DF_FILE

from sklearn.decomposition import PCA
# compute PCA on df_joined
from sklearn.preprocessing import StandardScaler

dict_gene_functions = {
            'polysaccharide': ['epsA', 'epsB', 'epsC', 'epsD', 'epsE', 'epsF', 'epsG',
                               'epsH', 'epsI', 'epsJ', 'epsK', 'epsL', 'epsM', 'epsN', 'epsO', 'galE'],
            'amyloid': ['tapA', 'sipW', 'tasA'],
            'regulation': ['SlrR', 'SlrA', 'SinR', 'SinI', 'KinC', 'KinD', 'Spo0A', 'PtkA', 'TkmA', 'PtpZ',
                           'DegU', 'DegQ', 'YmdB', 'FtsH', 'Veg', 'MstX', 'YugO'],
            'biofilm formation': ['AmpS', 'FloT', 'LuxS', 'RemA', 'RemB', 'Rny', 'Sfp/1', 'Sfp/2',
                          'SpeA', 'SpeD', 'SwrAA', 'YisP', 'YlbF', 'YmcA', 'YvcA', 'YwcC', 'YxaB'],
            'pellicle': ['Hag', 'FlgE', 'FliF', 'MotA', 'SigD', 'CheA',
                         'CheY', 'CheD', 'CheV', 'HemAT'],
            'secreted': ['abn2', 'abnA', 'amyE', 'aprE', 'artP', 'bdbA', 'bglC', 'blyA', 'bpr', 'bslA', 'csn',
                         'cwlO', 'cwlQ', 'dacA', 'dacC', 'dacF', 'efeB', 'epeX', 'epr', 'fecC', 'feuA', 'fhuD',
                         'flgB', 'flgC', 'flgE', 'flgK', 'flgM', 'flhO', 'flhP', 'fliD', 'fpbQ', 'ganB', 'ggt',
                         'glpQ', 'gmuG', 'hag', 'lip', 'lipB', 'ltaS', 'lytD', 'melE', 'metQ', 'mpr', 'nprB',
                         'nprE', 'oppA', 'pbpA', 'pbpB', 'pbpC', 'pbpX', 'pel', 'pelB', 'pelC', 'penP', 'pgdS',
                         'phoA', 'phoB', 'phoD', 'phy', 'rbsB', 'sacB', 'sacC', 'sivA', 'sivC', 'sunA', 'tagT',
                         'tasA', 'vpr', 'wapA', 'wprA', 'xepA', 'xkdG', 'xkdK', 'xkdM', 'xlyA', 'xlyB', 'xynA',
                         'xynC', 'xynD', 'yacD', 'ybdN', 'ybdO', 'ybfO', 'ybxI', 'ydaJ', 'ydhF', 'ydjM', 'ydjN',
                         'yfiY', 'yfkN', 'yhcR', 'yjcM', 'yjfA', 'yoaJ', 'yoaW', 'yocH', 'yolA', 'yolB', 'yqiH',
                         'yqxI', 'yrpD', 'yurI', 'yvgO', 'ywaD', 'ywoF', 'yxaL', 'yxeB', 'yxkC', 'zinT']
        }

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


def plot_all_features_heatmap(df, species_name, show=False):
    # plot heatmap of correlations
    df = df.replace(to_replace='None', value=np.nan).dropna()
    table = df.corr()
    fig = plt.figure()
    fig.set_size_inches(30, 30, forward=True)
    plt.rc('font', size=15)
    sns.heatmap(table, annot=False)

    if show:
        plt.show()
    out_file = '../../data/data_graphs/features_correlations/{}.png'.format(species_name)
    create_dir_if_not_exists(out_file)
    fig.savefig(out_file)


def plot_pca(df_labeled, target_column, sp_name, features_of_interest=[]):
    features_of_interest = ['NONPOLAR_AA',
                            'AROMATIC_AA',
                            'H1',
                            'MASS',
                            'P2',
                            'POSITIVE_AA',
                            'SASA']
    if len(features_of_interest) == 0:
        features = list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys())
    else:
        features = features_of_interest

    # Separating out the features
    df_labeled = df_labeled.dropna(0)  # NOTE: in BS3610 we had 58 genes of NAN in some of the features

    x = df_labeled.loc[:, features].values

    # Separating out the target
    y = df_labeled.loc[:, [target_column]].values

    # Standardizing the features
    x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x, )
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])

    print(len(principalDf))
    print(len(df_labeled))
    finalDf = pd.concat([principalDf, df_labeled[[target_column]]], axis=1)
    print(len(finalDf))
    finalDf = finalDf.dropna(0)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC1', fontsize=16)
    ax.set_ylabel('PC2', fontsize=16)
    ax.set_title(sp_name, fontsize=18)
    # targets = [False, True]
    # targets = list(set(finalDf[target_column].values))

    # remove the most frequent labels and add to the end
    from collections import Counter
    label_counter_sorted = Counter(finalDf[target_column].values).most_common()
    targets = [x[0] for x in label_counter_sorted]
    print(Counter(finalDf[target_column].values))

    colors = ['#bbdefb', '#5c6bc0','#ab47bc','#ef5350','#ffa726','#ffca28','#d4e157']

    # assert(len(targets) == len(colors))
    for i, (target, color) in enumerate(zip(targets, colors)):
        alpha = 0.3 if i == 0 else 0.8
        indicesToKeep = finalDf[target_column] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'],
                   finalDf.loc[indicesToKeep, 'principal component 2'],
                   s=50, c=color, alpha=alpha)
    # ax.legend(targets, fontsize=15, loc=1)
    ax.tick_params(axis='both', which='major', labelsize=13)
    right_side = ax.spines["right"]
    right_side.set_visible(False)

    top_side = ax.spines["top"]
    top_side.set_visible(False)

    out_file = '../../data/data_graphs/poster_pcas/{}.png'.format(sp_name)
    create_dir_if_not_exists(out_file)
    if sp_name == 'bacillus_subtilis':
        plt.title('{}\n({:.0f}% explained variance)'.format("Bacillus subtilis", 100*sum(pca.explained_variance_ratio_)), fontsize=25)
    elif sp_name == 'mega_species':
        plt.title('{}\n({:.0f}% explained variance)'.format("'Mega' Species", 100 * sum(pca.explained_variance_ratio_)), fontsize=25)
    else:
        plt.title('{}\n({:.0f}% explained variance)'.format(sp_name, 100 * sum(pca.explained_variance_ratio_)), fontsize=25)
    plt.grid(False)
    fig.savefig(out_file)


    # another trial with seaborn
    import seaborn as sns
    sns.set(style='white')

    # Set your custom color palette
    # sns.set_palette(sns.color_palette('pastel'))
    sns.set_palette(sns.color_palette(colors))

    sns.lmplot(x="principal component 1", y="principal component 2",
               data=finalDf,
               fit_reg=False,
               hue='is_gene_of_interest',  # color by cluster
               legend=True,
               scatter_kws={"s": 50, 'alpha': 0.5})  # specify the point size
    plt.show()

    # g = sns.scatterplot(x="principal component 1", y="principal component 2", hue="is_gene_of_interest",
    #                     data=finalDf, palette='colorblind',
    #                     legend='full')
    # # g.set(xscale="log")
    # plt.show()

# plot for each pair of species subplot of all features
def plot_histograms_pairwise_species():
    # create all pairs of species
    rows = 1
    cols = 10
    feature_idx = 0
    is_done = False

    #all_features = list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys())
    all_features = [
        'GC_CONTENT',
        'HYDROPHOBIC_AA',
        'NEGATIVE_AA',
        'NONPOLAR_AA',
        'MASS',
        'SASA',
        'P2',
        'V',
        'PI',
        'H1']
    for i in range(len(species_names)):
        species1 = species_names[i]
        goi_list = [item for sublist in list(dict_gene_functions.values()) for item in sublist]

        df1 = pd.read_pickle('../../data/data_outputs/features_' + species1 + '.pickle')
        df1['GENE_NAME'] = df1['GENE_NAME'].apply(lambda x: x.lower())
        df1 = df1[(df1['PRODUCT_TYPE'] == 'CDS') & (df1['IS_PSEUDO'] == False)]
        df1 = mark_genes_of_interest(df1, goi_list)
        df1 = df1[df1['is_gene_of_interest'] == True]

        print('df1 len: ', len(df1))
        for j in range(i + 1, len(species_names)):
            # plt.figure(figsize=(200, 200))
            is_done = False
            feature_idx = 0
            species2 = species_names[j]

            df2 = pd.read_pickle('../../data/data_outputs/features_' + species2 + '.pickle')
            df2['GENE_NAME'] = df2['GENE_NAME'].apply(lambda x: x.lower())
            df2 = df2[(df2['PRODUCT_TYPE'] == 'CDS') & (df2['IS_PSEUDO'] == False)]
            df2 = mark_genes_of_interest(df2, goi_list)
            df2 = df2[df2['is_gene_of_interest'] == True]
            print('df2 len: ', len(df2))
            print('Currently analyzing: ', species1, species2)
            # for each feature
            # for gene_feature in list(GeneFeatures.__members__.keys()) + list(ProteinFeatures.__members__.keys()):

            fig, axs = plt.subplots(rows, cols)
            fig.set_size_inches(60, 4, forward=True)

            plt.rc('font', size=20)
            fig.suptitle('{}({})-{}({})'.format(species1, len(df1), species2, len(df2)))

            for r in range(rows):
                if is_done:
                    break

                for c in range(cols):
                    gene_feature = all_features[feature_idx]
                    print('feature=', gene_feature, ' feature_idx=', feature_idx, ' out of:', len(all_features) - 1)
                    # axs[r, c].plot(x, y)
                    col1 = df1[gene_feature]
                    col2 = df2[gene_feature]

                    if rows == 1:
                        # specific code for plotting all of the subplots horizontally
                        col1[~col1.isnull()].plot.hist(bins=100, ax=axs[c], alpha=0.5, density=True)
                        col2[~col2.isnull()].plot.hist(bins=100, ax=axs[c], color='red', alpha=0.5, density=True)
                        axs[c].set_title(gene_feature)
                    else:
                        col1[~col1.isnull()].plot.hist(bins=100, ax=axs[r, c], alpha=0.5, density=True)
                        col2[~col2.isnull()].plot.hist(bins=100, ax=axs[r, c], color='red', alpha=0.5, density=True)
                        axs[r, c].set_title(gene_feature)

                    feature_idx += 1
                    if feature_idx == len(all_features):
                        is_done = True
                        break
            for ax in axs.flat:
                ax.set_ylabel('Counts', fontsize=22)
                # ax.set_xlabel('Values', fontsize=22)
                ax.tick_params(axis='both', which='major', labelsize=22)
                # ax.set(autoscale_on=True)

            # Hide x labels and tick labels for top plots and y ticks for right plots.
            for ax in axs.flat:
                ax.label_outer()
            # plt.show()
            out_file = '../../data/data_graphs/poster_histograms/{}_{}.png'.format(species1, species2)
            create_dir_if_not_exists(out_file)
            fig.savefig(out_file, dpi=300)
            # plt.show()



#
#     # plt.figure
#     # for each file in  pickle dir
#     # plt.subplot()

# #FOR POSTER GRAPHS ONLY
def plot_pca_one_species(df_species, species_name, add_function_labels=False, features_of_interest=[]):
    df_cds = df_species[df_species['PRODUCT_TYPE'] == 'CDS'].copy()

    # change to lower in order to join with the uniprot table
    df_cds['GENE_NAME'] = df_cds['GENE_NAME'].apply(lambda x: x.lower())

    use_uniprot = False
    if use_uniprot:
        # UNIPROT DF:
        UNIPROT_FILE = '../../data/data_inputs/bs168_uniprot.tab'
        df_UP = pd.read_csv(UNIPROT_FILE, sep='\t', header=0)
        df_UP = df_UP.fillna('')
        df_UP['isSecreted'] = df_UP['Subcellular location [CC]'].apply(lambda x: re.search('Secreted', x) is not None)
        df_UP['isExtracellular'] = df_UP['Topological domain'].apply(
            lambda x: re.search('Extracellular', x) is not None)
        # Note: we fetched only the 1st name in the list
        df_UP['GENE_NAME'] = df_UP['Gene names'].apply(
            lambda x: x.lower().split()[0] if x != '' else 'NA_GENE_NAME_UNIPROT')

        df_joined = pd.merge(df_cds, df_UP, on="GENE_NAME")

        df_joined['is_gene_of_interest'] = df_joined['isSecreted'] | df_joined['isExtracellular']
        df_genes_interest = df_joined[df_joined['isSecreted'] | df_joined['isExtracellular']]
        df_other_genes = df_joined[~df_joined['isSecreted'] & ~df_joined['isExtracellular']]
    else:
        df_joined = df_cds
        df_joined = df_joined.reset_index()

    # plot_all_features_histograms(df_genes_interest, suffix='_Secreted_Extracellular')
    # plot_all_features_histograms(df_other_genes, suffix='_Other')
    # plot_all_features_heatmap(df_genes_interest, suffix='_Secreted_Extracellular')
    # plot_all_features_heatmap(df_other_genes, suffix='_Other')

    # set all labels to true just in order to have a single color
    if not add_function_labels:
        df_joined['is_gene_of_interest'] = True  # for one coloring only

    # adding new labels from:
    # http://subtiwiki.uni-goettingen.de/wiki/index.php/Biofilm_formation#Key_genes_and_operons_involved_in_biofilm_formation
    else:
        df_joined['is_gene_of_interest'] = 'other'


        # epsA - epsB - epsC - epsD - epsE - epsF - epsG - epsH - epsI - epsJ - epsK - epsL - epsM - epsN - epsO
        # genes_biofilm = ['galE', 'tapA', 'sipW', 'tasA', 'bslA']
        dict_biofilm_lower = {}
        for k in dict_gene_functions.keys():
            dict_biofilm_lower[k] = [x.lower() for x in dict_gene_functions[k]]
        # genes_biofilm = [x.lower() for x in genes_biofilm]
        # print('len(genes_biofilm):', len(genes_biofilm))

        for k_group in dict_biofilm_lower.keys():
            count_found = 0
            for gene in dict_biofilm_lower[k_group]:
                count_found += sum(df_joined['GENE_NAME'] == gene)
                df_joined.loc[df_joined['GENE_NAME'] == gene, ['is_gene_of_interest']] = k_group
            print('for {}, found {}'.format(k_group, count_found))
        # print(sum(df_joined['is_gene_of_biofilm']))
    plot_pca(df_joined, target_column='is_gene_of_interest', sp_name=species_name,
             features_of_interest=features_of_interest)
    # 'is_gene_of_biofilm')
    # sres = set(s1_lower).intersection(set(s2_lower))


def mark_genes_of_interest(df, genes_of_interest):
    df['is_gene_of_interest'] = False
    for gene in genes_of_interest:
        df.loc[df['GENE_NAME'] == gene.lower(), ['is_gene_of_interest']] = True
    return df


# # ORIGINAL MAIN:
# if __name__ == '__main__':
#    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
#    df_all = pd.read_pickle(FEATURES_DF_FILE)
#    df_cds = df_all[df_all['PRODUCT_TYPE'] == 'CDS']
#
#    # create features histograms and heatmap
#    plot_all_features_histograms(df_cds)
#    plot_all_features_heatmap(df_cds)
#
#    # scatter plot:
#    #sns.jointplot(x='MELTING_POINT', y='GC_CONTENT', data=df_cds)
#
#    print('Done main_analysis')


# last main:
if __name__ == '__main__':
    # plot the features for all pairs of species - PART A of poster
    #plot_histograms_pairwise_species()
    #exit(0)

    df_mega = pd.DataFrame()
    to_plot_pca = True
    # plotting PCA for all species
    for i in range(len(species_names)):
        species_name = species_names[i]
        print('working on: {}'.format(species_name))
        FEATURES_DF_FILE = '../../data/data_outputs/features_{}.pickle'.format(species_name)
        df_species = pd.read_pickle(FEATURES_DF_FILE)

        plot_all_features_heatmap(df_species, species_name)

        if to_plot_pca:
            plot_pca_one_species(df_species=df_species, species_name=species_names[i], add_function_labels=True,
                                features_of_interest=[])

            # concatenate the dataframes to create a "mega-species"
            # Stack the DataFrames on top of each other
            df_mega = pd.concat([df_mega, df_species], axis=0)

    #print(len(df_mega))
    if to_plot_pca:
        plot_pca_one_species(df_species=df_mega, species_name='mega_species', add_function_labels=True,
                            features_of_interest=[])
