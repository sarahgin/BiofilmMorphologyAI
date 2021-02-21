import pickle
import matplotlib.pyplot as plt
NEXT_4_NT_FILENAME = 'next_4_nt.pickle'
if __name__ == '__main__':
    with open(NEXT_4_NT_FILENAME, 'rb') as handle:
        next_nt = pickle.load(handle)

    kmers = list(next_nt.keys())
    max_percents = []
    max_labels = []
    for k in next_nt.keys():
        current_kmers_next_nt_dict = next_nt[k]
        m = max(current_kmers_next_nt_dict.values())
        s = sum(current_kmers_next_nt_dict.values())
        max_percents.append(m/s)
        max_key = max(current_kmers_next_nt_dict, key=current_kmers_next_nt_dict.get)
        max_labels.append(max_key)

    z = range(len(max_percents))
    y = max_percents

    fig, ax = plt.subplots()
    ax.scatter(z, y)
    for i, txt in enumerate(max_labels):
        ax.annotate(kmers[i] + '-' + txt, (z[i], y[i]))

    plt.show()
    print('done')

    # print('Loading pickle file: {}...'.format(NEXT_4_NT_FILENAME))
    # next_nt_df = pd.read_pickle(NEXT_4_NT_FILENAME)
    # next_nt_df["SUM"] = next_nt_df[["A", "C", "G", "T"]].sum(axis=1)
    # next_nt_df["A_PER"] = next_nt_df["A"] / next_nt_df["SUM"] * 100
    # next_nt_df["C_PER"] = next_nt_df["C"] / next_nt_df["SUM"] * 100
    # next_nt_df["G_PER"] = next_nt_df["G"] / next_nt_df["SUM"] * 100
    # next_nt_df["T_PER"] = next_nt_df["T"] / next_nt_df["SUM"] * 100
    #
    # next_nt_df["MAX_PER"] = next_nt_df[["A_PER", "C_PER", "G_PER", "T_PER"]].max(axis=1)
    # plt.scatter(range(len(next_nt_df["MAX_PER"])), next_nt_df["MAX_PER"])
    # plt.show()
    # exit(0)

    def create_next_nucleotide_df(species_df, product_type: str, min_gene_length=0, max_gene_length=np.inf):
        # Step 1: create a dict with hexamers as keys and values
        # are lists of all normalized positions
        next_nucleotide_dict = {}
        count = 1

        # filter out step
        species_df = species_df[species_df['PRODUCT_TYPE'] == product_type]
        species_df = species_df[species_df['DNA_LENGTH'] >= min_gene_length]
        species_df = species_df[species_df['DNA_LENGTH'] <= max_gene_length]

        for index, row in species_df.iterrows():  # for each Gene
            print(count)
            count += 1

            gene_next_nucleotide_dict = row[GeneralFeatures.HEXAMER_NEXT_NUCLEOTIDE.name]
            for k in gene_next_nucleotide_dict:
                hex_next_nucleotide_dict = gene_next_nucleotide_dict[k]
                if k in next_nucleotide_dict:
                    next_nucleotide_dict[k] = merge_add_dicts(next_nucleotide_dict[k], hex_next_nucleotide_dict)
                else:
                    next_nucleotide_dict[k] = hex_next_nucleotide_dict
                # if k in next_nucleotide_dict:
                #    next_nucleotide_dict[k] = next_nucleotide_dict[k] + np.array(list(hex_next_nucleotide_dict.values()))
                # else:
                #    next_nucleotide_dict[k] = np.array(list(hex_next_nucleotide_dict.values()))
        # TODO fix next line!
        next_nucleotide_df = pd.DataFrame.from_dict(next_nucleotide_dict, orient='index', columns=['A', 'C', 'G', 'T'])

        return next_nucleotide_df
