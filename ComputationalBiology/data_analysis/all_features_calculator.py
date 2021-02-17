# from pytictoc import TicToc
# t = TicToc() #create instance of class
from enum import Enum

import pandas as pd

from ComputationalBiology.bio_general import Species
from ComputationalBiology.data_analysis.gene_features_calculator import *
from ComputationalBiology.data_analysis.general_features_calculator import *
from ComputationalBiology.data_analysis.protein_features_calculator import *


class GeneralFeatures(Enum):
    GENE_ID = 1
    GENE_NAME = 2
    TYPE = 3
    PRODUCT_TYPE = 4
    STRAND = 5
    PRODUCT_DESCRIPTION = 6
    HEXAMER_DICT = 7
    HEXAMER_NEXT_NUCLEOTIDE = 8
    # CODON_DICT = 9
    # TODO: add tss?


# features of the DNA sequence
class GeneFeatures(Enum):
    GC_CONTENT = 1
    DNA_LENGTH = 2


# features of the protein sequence
class ProteinFeatures(Enum):
    HYDROPHOBIC_AA = 1
    HYDROPHILIC_AA = 2
    POLAR_AA = 3
    AROMATIC_AA = 4
    POSITIVE_AA = 5
    NEGATIVE_AA = 6
    NONPOLAR_AA = 7
    AA_LENGTH = 8

class KmerFeatures(Enum):
    PREFIX_SUFFIX_DICT = 1

general_features_map = {
    GeneralFeatures.GENE_ID: get_gene_id,
    GeneralFeatures.GENE_NAME: get_gene_name,
    GeneralFeatures.TYPE: get_type,
    GeneralFeatures.PRODUCT_TYPE: get_product_type,
    GeneralFeatures.STRAND: get_strand,
    GeneralFeatures.PRODUCT_DESCRIPTION: get_product_description,
}

gene_features_map = {
    GeneFeatures.GC_CONTENT: compute_gc_content,
    GeneFeatures.DNA_LENGTH: compute_gene_length
}

protein_features_map = {
    ProteinFeatures.HYDROPHOBIC_AA: compute_hydrophobic_aa,
    ProteinFeatures.HYDROPHILIC_AA: compute_hydrophilic_aa,
    ProteinFeatures.POLAR_AA: compute_polar_aa,
    ProteinFeatures.AROMATIC_AA: compute_aromatic_aa,
    ProteinFeatures.POSITIVE_AA: compute_positive_aa,
    ProteinFeatures.NEGATIVE_AA: compute_negative_aa,
    ProteinFeatures.NONPOLAR_AA: compute_nonpolar_aa,
    ProteinFeatures.AA_LENGTH: compute_protein_length
}

kmer_features_map = {
    KmerFeatures.PREFIX_SUFFIX_DICT: compute_all_suffixes_given_prefix,

    #GeneralFeatures.HEXAMER_DICT: compute_hexamer_positions,
    #GeneralFeatures.HEXAMER_NEXT_NUCLEOTIDE: compute_hexamer_next_nucleotide,
    #GeneralFeatures.CODON_DICT: compute_codon_counts #TODO codon bias
}


def create_species_df(spp: Species):
    # TODO: when adding new features, no need to recompute
    """
    Go over all of the genes of the given species
    :param spp:
    :return:
    """
    # gene_type is a list containing the required type
    # if it is an empty list it means that we should not filter out any type

    list_of_dict = []
    for gene_key in spp.all_genes:
        gene = spp.all_genes[gene_key]

        #  dictionary fro current gene:
        current_features_dict = create_gene_features_dict(gene)
        list_of_dict.append(current_features_dict)

        if len(list_of_dict) % 100 == 0:
            print("gene num:", len(list_of_dict))

    df = pd.DataFrame(list_of_dict)
    df = df.set_index(GeneralFeatures.GENE_ID.name)

    return df


def create_gene_features_dict(gene: Gene):
    features_dict = {}
    for key in GeneralFeatures:
        func = general_features_map[key]
        features_dict[key.name] = func(gene)

    for key in GeneFeatures:
        func = gene_features_map[key]
        features_dict[key.name] = func(gene.coding_sequence)

    if gene.gene_product is not None:
        for key in ProteinFeatures:
            func = protein_features_map[key]
            features_dict[key.name] = func(gene.gene_product.translation)

    return features_dict


def create_gene_features_dict(gene: Gene, prefix_length_min,
                                          prefix_length_max,
                                          suffix_length_min,
                                          suffix_length_max):
    seq = gene.coding_sequence
    result = {}
    #we need to save the total number of times the prefix
    #appears so that we can divide the counts by this number to get percentages
    for pos in range(gene.codon_start, len(seq) - 1): #TODO handle positions
        print('pos' + str(pos))
        for curr_prefix_len in range(prefix_length_min, prefix_length_max + 1):
            curr_prefix = seq[pos: pos + curr_prefix_len]
            for curr_suffix_len in range(suffix_length_min, suffix_length_max + 1):
                if pos + curr_prefix_len + curr_suffix_len > len(seq):
                    continue
                curr_suffix = seq[pos + curr_prefix_len: pos + curr_prefix_len + curr_suffix_len]
                print('prefix:' + curr_prefix + ',suffix:' + curr_suffix)
                if curr_prefix not in result.keys():
                    result[curr_prefix] = {}
                if curr_suffix not in result[curr_prefix].keys():
                    result[curr_prefix][curr_suffix] = 1
                else:
                    result[curr_prefix][curr_suffix] += 1
    return result

if __name__ == '__main__':
    g = Gene(coding_sequence = 'AGTCGCCAATTT', start=0, end=12, strand=1, gene_type='CDS', name="dummy", qualifiers = None)
    res = create_gene_features_dict(g, 3, 5, 2, 4)
    print(res)