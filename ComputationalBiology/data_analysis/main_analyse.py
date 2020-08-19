import os

from BiofilmMorphologyAI.ComputationalBiology.fileutils.gene_bank_parser import read_genbank_file
from BiofilmMorphologyAI.ComputationalBiology.genetics.Gene import Gene

if __name__ == '__main__':
    # gene_bank_file = '../../data/GeneticData/BS3610.gb'
    # assert (os.path.exists(gene_bank_file))  # making sure that the path is valid
    #
    # record_gb = read_genbank_file(gene_bank_file)
    g = Gene(start=0, end=99, strand=1, sequence='', gene_product=None, name="", id="")
    print(g.is_coding())