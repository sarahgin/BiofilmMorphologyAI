"""
" The goal: to parse a GeneBank file and retrieve the information as Python objects
"""
import os
from Bio import SeqIO

from ComputationalBiology.biology_utils.Species import Species


def read_genbank_file(gene_bank_file: str):
    """
    Returns the 1st record in the given genebank file
    """
    with open(gene_bank_file, 'r') as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
        return record_gb


# def init_all_genes(record_gb):
#     all_genes = {}
#     features = record_gb.features
#     for f in features:
#         if f.type == 'CDS':
#             current_gene = Gene(start=,end=, strand=, name=, locus_tag=)
#     # start: int, end: int, strand, name: str, locus_tag: str)
#
#     return all_genes


# def init_species(record_gb):
#     all_genes = init_all_genes(record_gb.features)
#
#     #  create a new Species instance
#     species_obj = Species(species_name="", species_id="", all_genes={})
#     return species_obj


def get_stats(record_gb):
    type_stats = {}
    features = record_gb.features
    for i in range(1, len(features)):
        if features[i].type in type_stats:
            type_stats[features[i].type] += 1
        else:
            type_stats[features[i].type] = 1

    return type_stats


if __name__ == '__main__':
    genbank_file = '../../data/GeneticData/BS3610.gb'
    assert(os.path.exists(genbank_file))  # making sure that the path is valid
    record_gb = read_genbank_file(genbank_file)

    print(get_stats(record_gb))
    # {'gene': 4417, 'CDS': 4296, 'rRNA': 30, 'tRNA': 86, 'ncRNA': 4, 'regulatory': 26, 'tmRNA': 1}

    # convert the gene bank record to our objects:
    print(record_gb)

# TODO: operons?
# TODO: assert that locus_tag is unique
#
# # for cds: use codon_start, translation_table

