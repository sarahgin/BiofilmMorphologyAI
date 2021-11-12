print('hello')
import matplotlib.pyplot as plt
import pickle

# from ComputationalBiology.data_analysis.main_parser_features_calc import PREFIX_SUFFIX_DF_FILE

PREFIX_SUFFIX_DF_FILE = '../../data/data_outputs/prefix_suffix_dict_' + 'BS3610' + '.pickle'
print('Loading pickle file: {}...'.format(PREFIX_SUFFIX_DF_FILE))
with open(PREFIX_SUFFIX_DF_FILE, 'rb') as pickle_file:
    prefix_suffix_dict, prefix_count = pickle.load(file=pickle_file)

    # normalize
    for k in prefix_suffix_dict:
        for j in prefix_suffix_dict[k]:
            prefix_suffix_dict[k][j] /= prefix_count[k]
            prefix_suffix_dict[k][j] *= 100

    max_dict = {}
    for k in prefix_suffix_dict:
        d = prefix_suffix_dict[k]
        d_max_val = max(d.values())

        # taking list of car values in v
        v = list(d.values())
        # taking list of car keys in v
        vk = list(d.keys())
        d_max_str = (vk[v.index(max(v))])

        max_dict[k + '-' + d_max_str] = d_max_val

    #sorted_max_dict = {k: v for k1, v1 in sorted(max_dict.items(), key=lambda item: item[1])}

    z = list(range(len(max_dict.values())))
    y = list(max_dict.values())

    fig, ax = plt.subplots()
    ax.scatter(z, y)
    for i, txt in enumerate(max_dict.keys()):
        print(z[i])
        print(y[i])
        print(txt)
        ax.annotate(txt, (z[i], y[i]))

    plt.show()
exit(0)

#test_visualization_2
import matplotlib.pyplot as plt
import pickle

# from ComputationalBiology.data_analysis.main_parser_features_calc import PREFIX_SUFFIX_DF_FILE

#PREFIX_SUFFIX_DF_FILE = '../../data/data_outputs/prefix_suffix_dict_' + 'Bacillus-clausii' + '.pickle'


print('Loading pickle file: {}...'.format(PREFIX_SUFFIX_DF_FILE))
with open(PREFIX_SUFFIX_DF_FILE, 'rb') as pickle_file:
    prefix_suffix_dict, prefix_count = pickle.load(file=pickle_file)

    # normalize
    for k in prefix_suffix_dict:
        for j in prefix_suffix_dict[k]:
            prefix_suffix_dict[k][j] /= prefix_count[k]
            prefix_suffix_dict[k][j] *= 100
            if prefix_count[k] < 100:
                prefix_suffix_dict[k][j] = 0

    max_dict = {}
    for k in prefix_suffix_dict:
        d = prefix_suffix_dict[k]
        d_max_val = max(d.values())

        # taking list of car values in v
        v = list(d.values())
        # taking list of car keys in v
        vk = list(d.keys())
        d_max_str = (vk[v.index(max(v))])

        max_dict[k + '-' + d_max_str] = d_max_val

    #sorted_max_dict = {k: v for k1, v1 in sorted(max_dict.items(), key=lambda item: item[1])}

    z = list(range(len(max_dict.values())))
    y = list(max_dict.values())

    fig, ax = plt.subplots()
    ax.scatter(z, y)

    for i, txt in enumerate(max_dict.keys()):
        if y[i] > 0:
            ax.annotate(txt, (z[i], y[i]))

    plt.show()
exit(0)

print('Loading pickle file: {}...'.format(NEXT_NT_DF_FILE))
    next_nt_df = pd.read_pickle(NEXT_NT_DF_FILE)
    next_nt_df["SUM"] = next_nt_df[["A", "C", "G", "T"]].sum(axis=1)
    next_nt_df["A_PER"] = next_nt_df["A"] / next_nt_df["SUM"] * 100
    next_nt_df["C_PER"] = next_nt_df["C"] / next_nt_df["SUM"] * 100
    next_nt_df["G_PER"] = next_nt_df["G"] / next_nt_df["SUM"] * 100
    next_nt_df["T_PER"] = next_nt_df["T"] / next_nt_df["SUM"] * 100

    next_nt_df["MAX_PER"] = next_nt_df[["A_PER", "C_PER", "G_PER", "T_PER"]].max(axis=1)
    plt.scatter(range(len(next_nt_df["MAX_PER"])), next_nt_df["MAX_PER"])
    plt.show()
    exit(0)

    # Load data and prepare df_cds:
    print('Loading pickle file: {}...'.format(FEATURES_DF_FILE))
    df_all = pd.read_pickle(FEATURES_DF_FILE)
    df_cds = df_all[df_all['PRODUCT_TYPE'] == 'CDS']

    # Load kmer df file into df_kmers and create three new columns
    kmers_dict_list = pd.read_pickle(KMERS_DF_FILE)
    kmers_df = pd.DataFrame(kmers_dict_list)
    kmers_df['MEAN'] = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: np.mean(x))
    kmers_df['VAR'] = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: np.var(x))
    kmers_df['COUNT'] = kmers_df['RELATIVE_POSITIONS'].apply(lambda x: len(x))
    kmers_df['GC_CONTENT'] = kmers_df['KMER'].apply(lambda x: compute_gc_content(x))

    kmers_df['A'] = kmers_df['KMER'].apply(lambda x: compute_a_content(x))

    kmers_df['REL_POS_HIST_100'] = kmers_df['RELATIVE_POSITIONS'] \
        .apply(lambda x: np.histogram(x, bins=100))

    kmers_df['FIRST_BIN_BIAS'] = kmers_df['REL_POS_HIST_100'] \
        .apply(lambda x: np.divide(x[0][0], np.mean(x[0][1:])))

    kmers_df['LAST_BIN_BIAS'] = kmers_df['REL_POS_HIST_100'] \
        .apply(lambda x: np.divide(x[0][-1], np.mean(x[0][:-1])))

    # compute correlation between FIRST_BIN_BIAS and GC_CONTENT
    table = kmers_df[['COUNT', 'GC_CONTENT', 'FIRST_BIN_BIAS', 'LAST_BIN_BIAS', 'A']].corr()
    fig = plt.figure(figsize=(8.0, 5.0))
    sns.heatmap(table, annot=True)
    plt.show()

    # top frequent kmers position distribution box plot only
    kmers_df_sorted = kmers_df.sort_values(by='COUNT', ascending=False)
    plt.boxplot(kmers_df_sorted.iloc[range(-5, 5)]['RELATIVE_POSITIONS'])
    plt.show()

    # plot a histogram of all relative positions for most and least frequent kmers
    for i in range(10):
        most_frequent_pos_list = kmers_df_sorted.iloc[i]['RELATIVE_POSITIONS']
        plt.figure(i)
        plt.hist(most_frequent_pos_list, bins=100)
        plt.title(kmers_df_sorted.iloc[i]['KMER'] + ' #: ' + str(len(most_frequent_pos_list)))
        plt.show()

# boxplots of the three columns
# kmers_df.boxplot(column=['MEAN'])
# plt.title('MEAN')
# plt.show()
# kmers_df.boxplot(column=['VAR'])
# plt.title('VAR')
# plt.show()
# kmers_df.boxplot(column=['COUNT'])
# plt.title('COUNT')
# plt.show()
