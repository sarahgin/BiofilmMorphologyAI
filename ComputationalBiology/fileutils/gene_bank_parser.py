"""
" The goal: to parse a GeneBank file and retrieve the information as Python objects
"""
import os
from Bio import SeqIO

"""
" Returns the 1st record in the given genebank file
"""
def read_gene_bank_file(gene_bank_file):
    with open(gene_bank_file, 'r') as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
        return record_gb

# TODO:
# for f in features:
#     res = [(g1, p1), (g2, p2)...]


if __name__ == '__main__':
    gene_bank_file = '../../data/GeneticData/BS3610.gb'
    assert(os.path.exists(gene_bank_file))  # making sure that the path is valid
    record_gb = read_gene_bank_file(gene_bank_file)
    # print(record_gb.annotations)
    # print(type(record_gb.annotations))

    features = record_gb.features
    print(record_gb.features)
    print(type(record_gb.features))
