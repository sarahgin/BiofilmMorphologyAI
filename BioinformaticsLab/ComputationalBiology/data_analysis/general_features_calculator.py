from BioinformaticsLab.ComputationalBiology.bio_general.Gene import Gene
import re


def get_gene_id(gene: Gene):
    """

    :type gene: Gene object
    """
    return gene.get_id()


def get_gene_name(gene: Gene):
    return gene.name


def get_type(gene: Gene):
    return gene.type


def get_product_type(gene: Gene):
    return gene.gene_product.type if gene.gene_product is not None else None


def get_strand(gene: Gene):
    return gene.strand


def get_product_description(gene: Gene):
    return gene.gene_product.description if gene.gene_product is not None else None


def get_is_pseudo(gene: Gene):
    return gene.gene_product.is_pseudo if gene.gene_product is not None else None


def compute_gc_content(gene: Gene):
    """
    Given a genomic sequence (str), compute GC content.
    input: DNA sequence (str)
    output:%GC
    """
    dna_sequence = gene.coding_sequence
    return compute_gc_count(dna_sequence) / len(dna_sequence) * 100


def compute_gc_count(dna_sequence: str):
    matches = re.findall(r'[GC]', dna_sequence.upper())
    return len(matches)


def compute_gene_length(gene: Gene):
    return len(gene.coding_sequence)
