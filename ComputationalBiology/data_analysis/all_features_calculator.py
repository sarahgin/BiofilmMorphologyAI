# from pytictoc import TicToc
# t = TicToc() #create instance of class
from enum import Enum

import pandas as pd

from ComputationalBiology.bio_general import Species
from ComputationalBiology.bio_general.bio_macros import SUFFIX_LENGTH_MAX, PREFIX_LENGTH_MAX, PREFIX_LENGTH_MIN, \
    SUFFIX_LENGTH_MIN
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
    #HEXAMER_DICT = 7
    #HEXAMER_NEXT_NUCLEOTIDE = 8
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
    # chemical features
    H1 = 9
    H2 = 10
    H3 = 11
    V = 12
    P1 = 13
    P2 = 14
    SASA = 15
    NCI = 16
    MASS = 17


class KmerFeatures(Enum):
    PREFIX_SUFFIX_DICT = 1


features_to_compute = {}

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
    ProteinFeatures.AA_LENGTH: compute_protein_length,
    #chemical properties
    ProteinFeatures.H1: compute_H1,
    ProteinFeatures.H2: compute_H2,
    ProteinFeatures.H3: compute_H3,
    ProteinFeatures.V: compute_V,
    ProteinFeatures.P1: compute_P1,
    ProteinFeatures.P2: compute_P2,
    ProteinFeatures.SASA: compute_SASA,
    ProteinFeatures.NCI: compute_NCI,
    ProteinFeatures.MASS: compute_MASS
}

kmer_features_map = {
    KmerFeatures.PREFIX_SUFFIX_DICT: create_gene_prefix_suffix_dict,

    # GeneralFeatures.HEXAMER_DICT: compute_hexamer_positions,
    # GeneralFeatures.HEXAMER_NEXT_NUCLEOTIDE: compute_hexamer_next_nucleotide,
    # GeneralFeatures.CODON_DICT: compute_codon_counts #TODO codon bias
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

        if len(list_of_dict) % 1000 == 0:
            print("create_species_df: gene num-", len(list_of_dict))

    df = pd.DataFrame(list_of_dict)
    df = df.set_index(GeneralFeatures.GENE_ID.name)

    return df


def create_gene_features_dict(gene: Gene):
    features_dict = {}
    for key in GeneralFeatures:
        if key in features_to_compute or len(features_to_compute) == 0:
            func = general_features_map[key]
            features_dict[key.name] = func(gene)

    for key in GeneFeatures:
        if key in features_to_compute or len(features_to_compute) == 0:
            func = gene_features_map[key]
            features_dict[key.name] = func(gene.coding_sequence)

    #TEMPORARILY COMMENTED OUT
    #for key in KmerFeatures:
    #    if key in features_to_compute or len(features_to_compute) == 0:
    #        func = kmer_features_map[key]
    #        features_dict[key.name] = func(gene,
    #                                       PREFIX_LENGTH_MIN,
    #                                       PREFIX_LENGTH_MAX,
    #                                       SUFFIX_LENGTH_MIN,
    #                                       SUFFIX_LENGTH_MAX)

    if gene.gene_product is not None:
        for key in ProteinFeatures:
            if key in features_to_compute or len(features_to_compute) == 0:
                func = protein_features_map[key]
                features_dict[key.name] = func(gene.gene_product.translation)

    return features_dict

# if __name__ == '__main__':
#    g = Gene(coding_sequence = 'AGTCGCCAATTT', start=0, end=12, strand=1, gene_type='CDS', name="dummy", qualifiers = None)
#    res, counts = create_gene_features_dict(g, 3, 5, 2, 4)
#    print(res)
