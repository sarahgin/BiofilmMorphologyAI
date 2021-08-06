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
    print('e')
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    sequences.append(str(seq_record.seq))



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
print(clusterid)
print(np.unique(clusterid))

max_found = 0
max_id = 0
for i in range(0, 1000):
    num_curr = sum(clusterid==i)
    if num_curr > max_found:
        max_found = num_curr
        max_id = i
print(num_curr, max_id)
print('done')