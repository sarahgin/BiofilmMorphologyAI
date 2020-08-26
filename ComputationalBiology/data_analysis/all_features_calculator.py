from enum import Enum
import pandas as pd

from ComputationalBiology.biology_utils import Species
from ComputationalBiology.biology_utils.Gene import Gene
from ComputationalBiology.data_analysis.gene_features_calculator import *
from ComputationalBiology.data_analysis.protein_features_calculator import *


class GeneFeatures(Enum):
    GC_CONTENT = 1
    MELTING_POINT = 2
    LENGTH = 3


class ProteinFeatures(Enum):
    HYDROPHOBIC_AA = 1
    HYDROPHILIC_AA = 2
    POLAR_AA = 3
    AROMATIC_AA = 4
    POSITIVE_AA = 5
    NEGATIVE_AA = 6
    NONPOLAR_AA = 7
    LENGTH = 8


gene_features_map = {
        GeneFeatures.GC_CONTENT: compute_gc_content,
        GeneFeatures.MELTING_POINT: compute_melting_point,
        GeneFeatures.LENGTH: compute_gene_length
        }

protein_features_map = {
        ProteinFeatures.HYDROPHOBIC_AA: compute_hydrophobic_aa,
        ProteinFeatures.HYDROPHILIC_AA: compute_hydrophilic_aa,
        ProteinFeatures.POLAR_AA: compute_polar_aa,
        ProteinFeatures.AROMATIC_AA: compute_aromatic_aa,
        ProteinFeatures.POSITIVE_AA: compute_positive_aa,
        ProteinFeatures.NEGATIVE_AA: compute_negative_aa,
        ProteinFeatures.NONPOLAR_AA: compute_nonpolar_aa,
        ProteinFeatures.LENGTH: compute_protein_length
}


def calculate_all_features_species(spp: Species):
    # gene_type is a list containing the required type
    # if it is an empty list it means that we should not filter out any type
    df = pd.DataFrame()
    for gene_key in spp.all_genes:
        gene = spp.all_genes[gene_key]

        current_features_dict = calculate_all_features_gene(gene)
        df = df.append(current_features_dict, ignore_index=True)
        # break

    print(df)
    return df


def calculate_all_features_gene(gene: Gene):

    name = gene.qualifiers['gene'][0] if 'gene' in gene.qualifiers else ''
    product_type = gene.gene_product.type if gene.gene_product is not None else ''
    product_dsc = gene.gene_product.qualifiers['product'] if gene.gene_product is not None else ''

    features_dict = {'GENE_ID': gene.get_id(), 'GENE_NAME': name, 'TYPE': gene.type,
                     'PRODUCT_TYPE': product_type, 'STRAND': gene.strand,
                     'PRODUCT_DESCRIPTION': product_dsc}

    for gene_feature in GeneFeatures:
        key = gene_feature
        func = gene_features_map[key]
        features_dict[key.name] = func(gene.coding_sequence)
    return features_dict


