from Bio import SeqIO
import os
import numpy as np
from Bio.Cluster import kcluster
import pandas as pd


def compute_concensus(sub_list_seqs):
    return sub_list_seqs[0]


def report_clusters_to_file(sequences, clusterid, nclusters, fname):
    assert(len(sequences) == len(clusterid))
    data = []
    for x in range(0, nclusters):
        filter = clusterid == x
        filtered_list = [i for indx,i in enumerate(sequences) if filter[indx] == True]

        if len(filtered_list) == 1:
            data.append([x, 'Singelton', filtered_list[0]])
        else:
            peptide_num = 1
            data.append([x, 'Consnesus', compute_concensus(filtered_list)])
            for seq in filtered_list:
                data.append([x, peptide_num, seq])
                peptide_num += 1
    df = pd.DataFrame(data, columns=['Cluster Number', 'Peptide Number', 'Peptide'])
    df = df.astype(str)
    df.to_csv(fname, index=False, sep='\t')


if __name__ == '__main__':
    # Read Fasta file:
    organism = 'human1'
    dir_path = 'data//' + organism + '//'
    fasta_file = dir_path + 'top_subseqs_' + organism + '.txt'
    assert (os.path.exists(fasta_file))

    sequences = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(seq_record.seq))

    print('num of seqs: ', len(sequences))

    # sequences = sequences[:20]

    # Create clusters
    matrix = np.asarray([np.frombuffer(s.encode(),'u1') for s in sequences])
    nclusters = 1000
    npass = 500
    clusterid, error, nfound = kcluster(matrix, nclusters=nclusters, npass=npass, method='a')
    print('error: ', error)
    print('nfound: ', nfound)

    # Report result to file:
    fname = 'data/human1/noas_clusters_k_{}_total_{}_npass{}.csv'.format(nclusters,  len(sequences), npass)
    report_clusters_to_file(sequences, clusterid, nclusters, fname)

    print('done')
    exit()

