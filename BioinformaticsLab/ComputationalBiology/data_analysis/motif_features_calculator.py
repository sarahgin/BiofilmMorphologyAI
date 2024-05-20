import re
from BioinformaticsLab.ComputationalBiology.bio_general.Gene import Gene
import BioinformaticsLab.ComputationalBiology.biore.biore as biore
import numpy as np


def has_dna_motif(gene: Gene, motif: str):
    dna_sequence = gene.coding_sequence
    res = re.search(motif, dna_sequence, flags=re.IGNORECASE)
    return res is not None


def get_dna_motif_relative_positions(gene: Gene, motif: str):
    dna_sequence = gene.coding_sequence
    res = re.finditer(motif, dna_sequence, flags=re.IGNORECASE)
    start_positions = [np.round(e.span()[0] / len(dna_sequence), 2) for e in res]
    # TODO: if there are no matches at all, start_positions is an empty list
    # this will be a problem for np.mean(start_positions)
    return start_positions


def has_tataat_motif(gene: Gene):
    return has_dna_motif(gene, motif='tataat')


def get_tataat_relative_positions(gene: Gene):
    return get_dna_motif_relative_positions(gene, motif='tataat')


def has_tttatt_motif(gene: Gene):
    return has_dna_motif(gene, motif='tttatt')


def get_tttatt_relative_positions(gene: Gene):
    return get_dna_motif_relative_positions(gene, motif='tttatt')


#### Protein motifs ####
def has_protein_motif(gene: Gene, motif: str):
    '''
    motif is a biore motif (i.e. a letter represents a group,
    escaped letter is an aa) --> uses biore.search
    '''
    aa_sequence = gene.gene_product.translation
    res = biore.search(motif, aa_sequence, flags=re.IGNORECASE)
    return res is not None


def get_protein_motif_relative_positions(gene: Gene, motif: str):
    '''
    motif is a biore motif (i.e. a letter represents a group,
    escaped letter is an aa) --> uses biore.finditer
    '''
    aa_sequence = gene.gene_product.translation
    res = biore.finditer(motif, aa_sequence, flags=re.IGNORECASE)
    start_positions = [np.round(e.span()[0] / len(aa_sequence), 2) for e in res]
    # TODO: if there are no matches at all, start_positions is an empty list
    # this will be a problem for np.mean(start_positions)
    return start_positions


def has_ppkl_motif(gene: Gene):
    # TODO: this function is relevant only for gene with protein sequence
    return has_protein_motif(gene, '\P\P\P.{1,3}\L\L\L')


def get_ppkl_relative_positions(gene: Gene):
    return get_protein_motif_relative_positions(gene, '\P\P\P.{1,3}\L\L\L')
