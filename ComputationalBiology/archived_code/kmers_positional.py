# #filter out all positions that are greater than 100 nucleotides from gene beginning
#     kmers_df['ABSOLUTE_POSITIONS_FIRST_100'] \
#         = kmers_df['ABSOLUTE_POSITIONS'].apply(lambda x: [y for y in x if y <= 100])
#     #convert into 10-length bin vectors for histogram
#     kmers_df['POSITIONAL_BINS_DISTRIBUTION'] \
#         = kmers_df['ABSOLUTE_POSITIONS_FIRST_100'].apply(lambda x: np.histogram(x, bins=10)[0])
#     #save all kmer histogram pngs (4096)
#     kmers_pngs = kmers_df['ABSOLUTE_POSITIONS_FIRST_100'].to_list()
#     for index, row in kmers_df.iterrows():
#         kmer_name = row['KMER']
#         values = row['POSITIONAL_BINS_DISTRIBUTION']
#         fig = plt.figure()
#         plt.bar(range(10), values)
#         plt.title(kmer_name)
#         fig.savefig('../../data/data_outputs/PNGS/' + kmer_name + '.png')
#         plt.close()
#     exit(0)