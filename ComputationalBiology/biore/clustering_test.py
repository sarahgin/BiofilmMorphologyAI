from Bio.Cluster import kcluster
import numpy as np
from Bio import SeqIO
import os
# read seqs from fasta
organism = 'human1'
dir_path = 'data//' + organism + '//'

fasta_file = dir_path + 'top_subseqs_' + organism + '.txt'

assert (os.path.exists(fasta_file))
sequences = []


from Bio import SeqIO
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    # print('e')
    # print(seq_record.id)
    # print(repr(seq_record.seq))
    # print(len(seq_record))
    sequences.append(str(seq_record.seq))

print('num of seqs: ', len(sequences))

# sequences = [
# 'ADHAMKCAIROSURBANDJVUGLOBALIZATIONANDURBANFANTASIESPLA',
#      'ADHAMKCAIROSURBANDJVUGLOBALIZATIONANDURBANFANTASIESPLA',
#      'AGGESTAMKTHEARABSTATEANDNEOLIBERALGLOBALIZATIONTHEARAB',
#      'AGGESTAMKTHEARABSTATEANDNEOLIBERALGLOBALIZATIONTHEARAB',
#      'AGGESTAMKTHEARABSTATEANDNEOLIBERALGLOBALIZATIONTHEARAB'
#  ]
matrix = np.asarray([np.frombuffer(s.encode(),'u1') for s in sequences])
# matrix = np.asarray([np.fromstring(s, dtype=np.uint8) for s in sequences])
clusterid, error, nfound = kcluster(matrix, nclusters=1000)
# print(clusterid)
# print(np.unique(clusterid))

max_found = 0
max_id = 0
for x in range(0, 1000):
    num_curr = sum(clusterid==x)
    if num_curr > max_found:
        max_found = num_curr
        max_id = x
print(max_found, max_id)

print('Printing top cluster:')
filter = clusterid == max_id
filtered_list = [i for indx,i in enumerate(sequences) if filter[indx] == True]
for s in filtered_list:
    print(s)

n = 10
print('Printing all clusters with at least {} samples:'.format(n))
for x in range(0, 1000):
    num_curr = sum(clusterid==x)
    if num_curr < n:
        continue
    filter = clusterid == x
    filtered_list = [i for indx, i in enumerate(sequences) if filter[indx] == True]
    for s in filtered_list:
        print(s)
    print('--------------')
print('done')