"""
" The goal: to parse a GeneBank file and retrieve the information as Python objects
"""
import os
from Bio import SeqIO
# importing "collections" for namedtuple()
import collections

# defining a new type

from ComputationalBiology.biology_utils.Gene import Gene

genbank_item_tuple = collections.namedtuple('genbank_feature', ['gene', 'gene_product'])


def read_genbank_file(gene_bank_file):
    """
    Returns the 1st record in the given genebank file
    """
    with open(gene_bank_file, 'r') as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
        return record_gb


def filter_tuples(list_tuples, product_type='CDS'):
    """
    Given a list of tuples of (genbank_feature)s, get a list of only the tuples such that
    the gene product type is the given product_type
    """
    result = []
    for t in list_tuples:
        if t.gene_product.type == product_type:
            result.append(t)

    return result


def compare_genomic_position(seq_pos1, seq_pos2):
    """
    type: (SeqFeature, SeqFeature) -> bool
    """
    return (seq_pos1.strand == seq_pos1.strand and
            seq_pos1.location.start == seq_pos2.location.start and
            seq_pos1.location.end == seq_pos2.location.end)


def genebank_to_species(gene_bank_file):
    record_gb = read_genbank_file(gene_bank_file) # get 1st record

    dict_all_genes = {}
    features = record_gb.features
    assert(features[0].type == 'source')  # TODO: add check instead of assertion
    for i in range(1, len(features)):
        if features[i].type == 'gene':
            # get gene id
            # dict_all_genes[gene_id] = Gene()
            Gene(start=0, end=99, strand=1, sequence='', gene_product=None, name="", id="")


def genbank_to_tuples(record_gb):
    """
    Given a gene bank record containing information of both genes and, return a list of
    named tuples containing the CDS and the protein info of each matching gene.
    """
    result = []
    all_types = set()

    # fetch the features info
    features = record_gb.features
    assert(features[0].type == 'source')  # TODO: add check instead of assertion
    for i in range(1, len(features)-1):
        # assert(features[i].type == 'gene'), i # TODO: add check instead of assertion
        # assert (features[i+1].type == 'CDS'),  i  # TODO: add check instead of assertion

        if features[i].type == 'regulatory':
            print(i)
            # print('features[i-1]: ',  features[i-1])
            print('features[i]: ', features[i])
            print()

        if features[i].type == 'gene':
            start_s = features[i].location.start.position
            end_s = features[i].location.end.position

            start_next = features[i+1].location.start.position
            end_next = features[i+1].location.end.position

            if start_s != start_next or end_s != end_next:
                print(i)
                print('features[i]:', features[i])
                print('features[i+1]:', features[i+1])

            # assert(compare_genomic_position(features[i], features[i+1]))
            # if features[i+1].type == 'CDS':
                # assert(features[i].qualifiers['gene'] == features[i+1].qualifiers['gene']), i
            result.append(genbank_item_tuple(features[i], features[i + 1]))

        all_types.add(features[i].type)
        all_types.add(features[i+1].type)

    # assert(len(result) == (len(features) - 1)/2)
    return result, all_types


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

    type_stats = get_stats(record_gb)
    print(type_stats)
    # print(record_gb.annotations)
    # print(type(record_gb.annotations))

    # features = record_gb.features
    # print(record_gb.features)
    #
    lst, all_types = genbank_to_tuples(record_gb)
    # print(lst)
    # print(all_types)
    #
    # lst_filtered = filter_tuples(lst, product_type='CDS')
    # # TODO: compare gene name, position
    # print(lst_filtered)
    # for t in lst_filtered:
    #     print(t.gene)
    #     print(t.gene_product)
    #     print(' ')
