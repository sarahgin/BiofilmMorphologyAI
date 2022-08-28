# from pytictoc import TicToc
# t = TicToc() #create instance of class
from enum import Enum

from BioinformaticsLab.ComputationalBiology.bio_general import Species
from BioinformaticsLab.ComputationalBiology.data_analysis.discovery_calculator import *
from BioinformaticsLab.ComputationalBiology.data_analysis.general_features_calculator import *
from BioinformaticsLab.ComputationalBiology.data_analysis.motif_features_calculator import has_tataat_motif, \
    get_tataat_relative_positions, has_tttatt_motif, get_tttatt_relative_positions, has_ppkl_motif, \
    get_ppkl_relative_positions
from BioinformaticsLab.ComputationalBiology.data_analysis.protein_features_calculator import *


class GeneralFeatures(Enum):
    GENE_ID = 1
    GENE_NAME = 2
    TYPE = 3
    PRODUCT_TYPE = 4
    STRAND = 5
    PRODUCT_DESCRIPTION = 6
    IS_PSEUDO = 7
    GC_CONTENT = 8
    GENE_LENGTH = 9
    # TODO: add tss? (transcription start site)


# features of the protein sequence
class ProteinFeatures(Enum):
    HYDROPHOBIC_AA = 1
    HYDROPHILIC_AA = 2
    POLAR_AA = 3
    AROMATIC_AA = 4
    POSITIVE_AA = 5
    NEGATIVE_AA = 6
    NONPOLAR_AA = 7
    PROTEIN_LENGTH = 8

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

    #pH-related:
    PKA_COOH = 18
    PKA_NH = 19,
    PI = 20


class DNAMotifFeatures(Enum):
    HAS_TATAAT = 1
    RELATIVE_POSITIONS_TATAAT = 2
    HAS_TTTATT = 3
    RELATIVE_POSITIONS_TTTATT = 4


class ProteinMotifFeautres(Enum):
    HAS_PPKL = 1  # TODO: this is  a made up motif
    RELATIVE_POSITIONS_PPKL = 2  # TODO: this is  a made up motif


class PrefixSuffixFeatures(Enum):
    # represents the 3 dictionaries: full_dict, counts_dict
    # and relative_positions_dict
    GENE_PREFIX_SUFFIX_TUPLE = 1
    # PROTEIN_PREFIX_SUFFIX_TUPLE = 2


features_to_compute = {}


prefix_suffix_map = {
    PrefixSuffixFeatures.GENE_PREFIX_SUFFIX_TUPLE: create_gene_prefix_suffix_dict,
    # PrefixSuffix.PROTEIN_PREFIX_SUFFIX_TUPLE: create_protein_prefix_suffix_dict
}

general_features_map = {
    GeneralFeatures.GENE_ID: get_gene_id,
    GeneralFeatures.GENE_NAME: get_gene_name,
    GeneralFeatures.TYPE: get_type,
    GeneralFeatures.PRODUCT_TYPE: get_product_type,
    GeneralFeatures.STRAND: get_strand,
    GeneralFeatures.PRODUCT_DESCRIPTION: get_product_description,
    GeneralFeatures.IS_PSEUDO: get_is_pseudo,
    GeneralFeatures.GC_CONTENT: compute_gc_content,
    GeneralFeatures.GENE_LENGTH: compute_gene_length,
}


protein_features_map = {
    ProteinFeatures.HYDROPHOBIC_AA: compute_hydrophobic_aa,
    ProteinFeatures.HYDROPHILIC_AA: compute_hydrophilic_aa,
    ProteinFeatures.POLAR_AA: compute_polar_aa,
    ProteinFeatures.AROMATIC_AA: compute_aromatic_aa,
    ProteinFeatures.POSITIVE_AA: compute_positive_aa,
    ProteinFeatures.NEGATIVE_AA: compute_negative_aa,
    ProteinFeatures.NONPOLAR_AA: compute_nonpolar_aa,
    ProteinFeatures.PROTEIN_LENGTH: compute_protein_length,
    #Chemical properties:
    ProteinFeatures.H1: compute_H1,
    ProteinFeatures.H2: compute_H2,
    ProteinFeatures.H3: compute_H3,
    ProteinFeatures.V: compute_V,
    ProteinFeatures.P1: compute_P1,
    ProteinFeatures.P2: compute_P2,
    ProteinFeatures.SASA: compute_SASA,
    ProteinFeatures.NCI: compute_NCI,
    ProteinFeatures.MASS: compute_MASS,
    ProteinFeatures.PKA_COOH: compute_PKA_COOH,
    ProteinFeatures.PKA_NH: compute_PKA_NH,
    ProteinFeatures.PI: compute_PI
}

dna_motif_features_map = {
    DNAMotifFeatures.HAS_TATAAT: has_tataat_motif,
    DNAMotifFeatures.RELATIVE_POSITIONS_TATAAT: get_tataat_relative_positions,
    DNAMotifFeatures.HAS_TTTATT: has_tttatt_motif,
    DNAMotifFeatures.RELATIVE_POSITIONS_TTTATT: get_tttatt_relative_positions,
}

protein_motif_features_map = {
    ProteinMotifFeautres.HAS_PPKL: has_ppkl_motif, # TODO: this is  a made up motif
    ProteinMotifFeautres.RELATIVE_POSITIONS_PPKL: get_ppkl_relative_positions  # TODO: this is  a made up motif
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
    # threding from for to func

    list_of_dict = []
    for gene_key in spp.all_genes:
        gene = spp.all_genes[gene_key]

        #  dictionary for current gene:
        current_features_dict = create_gene_features_dict(gene)
        list_of_dict.append(current_features_dict)

        if len(list_of_dict) % 1000 == 0:
            print("create_species_df: gene num-", len(list_of_dict))

    df = pd.DataFrame(list_of_dict)
    df = df.set_index(GeneralFeatures.GENE_ID.name)

    return df


def create_gene_features_dict(gene: Gene):
    if gene.coding_sequence is None or gene.coding_sequence == '':
        print(gene)

    features_dict = {}
    # Computation for DNA sequences
    for key in GeneralFeatures:
        if key in features_to_compute or len(features_to_compute) == 0:
            func = general_features_map[key]
            features_dict[key.name] = func(gene)

    for key in DNAMotifFeatures:
        if key in features_to_compute or len(features_to_compute) == 0:
            func = dna_motif_features_map[key]
            features_dict[key.name] = func(gene)

    for key in PrefixSuffixFeatures:
        if key in features_to_compute or len(features_to_compute) == 0:
            func = prefix_suffix_map[key]
            features_dict[key.name] = func(gene)

    # Computation for protein sequences
    if gene.gene_product is not None:
        for key in ProteinFeatures:
            if key in features_to_compute or len(features_to_compute) == 0:
                func = protein_features_map[key]
                features_dict[key.name] = func(gene.gene_product.translation)

        for key in ProteinMotifFeautres:
            if key in features_to_compute or len(features_to_compute) == 0:
                func = protein_motif_features_map[key]
                features_dict[key.name] = func(gene)

    return features_dict

