from enum import Enum
import pandas as pd

from ComputationalBiology.biology_utils import Species
from ComputationalBiology.biology_utils.Gene import Gene
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


class GeneFeatures(Enum):
    GC_CONTENT = 1
    DNA_LENGTH = 2


class ProteinFeatures(Enum):
    HYDROPHOBIC_AA = 1
    HYDROPHILIC_AA = 2
    POLAR_AA = 3
    AROMATIC_AA = 4
    POSITIVE_AA = 5
    NEGATIVE_AA = 6
    NONPOLAR_AA = 7
    AA_LENGTH = 8


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

general_features_map = {
    GeneralFeatures.GENE_ID: get_gene_id,
    GeneralFeatures.GENE_NAME: get_gene_name,
    GeneralFeatures.TYPE: get_type,
    GeneralFeatures.PRODUCT_TYPE: get_product_type,
    GeneralFeatures.STRAND: get_strand,
    GeneralFeatures.PRODUCT_DESCRIPTION: get_product_description
}


def create_species_df(spp: Species):
    # TODO: when adding new features, no need to recompute
    # all of the table, just add the new feture
    """
    Go over all of the genes of the given species
    :param spp:
    :return:
    """
    # gene_type is a list containing the required type
    # if it is an empty list it means that we should not filter out any type
    df = pd.DataFrame()
    for gene_key in spp.all_genes:
        gene = spp.all_genes[gene_key]

        current_features_dict = create_gene_features_dict(gene)
        df = df.append(current_features_dict, ignore_index=True)

        if len(df) % 100 == 0:
            print("gene num:", len(df))

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
