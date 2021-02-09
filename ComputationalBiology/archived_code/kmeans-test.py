#   #TODO: move to machine learning file
#     #KMEANS: The following code is for clustering analysis only
#     array_to_cluster_list = kmers_df['POSITIONAL_BINS_DISTRIBUTION'].to_list()
#     #kmeans = KMeans(n_clusters=2, random_state=0).fit(array_to_cluster.to_list())
#
#     #Elbow method
#     distortions = []
#     silhouette_scores = []
#     K = range(2, 20)
#     for k in K:
#         print(k)
#         kmeanModel = KMeans(n_clusters=k, random_state=0).fit(array_to_cluster_list)
#
#         score = sklearn.metrics.silhouette_score(array_to_cluster_list, kmeanModel.labels_)
#
#         silhouette_scores.append(score)
#         distortions.append(sum(np.min(cdist(array_to_cluster_list, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) /
#                            len(array_to_cluster_list))
#
#     # Plot the elbow
#     plt.plot(K, distortions, 'bx-')
#     plt.plot(K, silhouette_scores, 'rx-')
#     plt.xlabel('k')
#     plt.ylabel('Distortion')
#     plt.title('The Elbow Method showing the optimal k')
#     plt.show()
#
#
#     print('done')
#
# # TODO: operons?
# # TODO: assert that locus_tag is unique
# #
# # # for cds: use codon_start, translation_table

#hex_norm = np.divide(hex_dict[k], row[GeneFeatures.DNA_LENGTH.name]-len(k)+1)