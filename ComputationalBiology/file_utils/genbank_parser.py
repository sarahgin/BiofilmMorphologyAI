"""
" The goal: to parse a GeneBank file and retrieve the information as Python objects
"""
import os
from Bio import SeqIO

from ComputationalBiology.biology_utils.Gene import Gene
from ComputationalBiology.biology_utils.GeneProduct import GeneProduct
from ComputationalBiology.biology_utils.Species import Species


def read_genbank_file(gene_bank_file: str):
    """
    Returns the 1st record in the given genebank file
    """
    with open(gene_bank_file, 'r') as input_handle:
        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
        return record_gb


def get_coding_gene_seq(genome_sequence, start, end, strand):

    assert(0 <= start < len(genome_sequence))
    assert(0 <= end < len(genome_sequence))

    seq = genome_sequence[start:end]
    if strand == -1:
        seq = seq.reverse_complement()

    return str(seq)


def init_all_genes(record_gb):
    features_list = record_gb.features
    genome_sequence = record_gb.seq

    all_genes = {}
    # for f in features_list:
    for i in range(len(features_list)) :
        f = features_list[i]
        #if f.type == 'CDS':

        if f.type == 'source':
            continue

        start = f.location.start
        end = f.location.end
        strand = f.location.strand
        coding_seq = get_coding_gene_seq(genome_sequence, start, end, strand)
        gene_key = str(start) + '_' + str(end) + '_' + str(strand)

        #ASSUMING EACH GENE AND PRODUCT HAVE THE SAME START+END+STRAND - TODO:check!
        if gene_key not in all_genes.keys():

            if 'gene' in f.qualifiers.keys():
                assert(len(f.qualifiers['gene']) <= 1), 'ERROR_' + features_list[i].type + \
                                                        '_' + features_list[i+1].type

            if f.type != 'gene' and f.type != 'regulatory':
                print('NO previous gene found for: ' + gene_key)
                continue

            name = f.qualifiers['gene'][0] if 'gene' in f.qualifiers else ''

            g = Gene(start=start,
                     end=end,
                     strand=strand,
                     gene_type=f.type,
                     name=name,
                     qualifiers=f.qualifiers,
                     coding_sequence=coding_seq)
            all_genes[gene_key] = g
        else:
            translation = '' if 'translation' not in f.qualifiers.keys() else f.qualifiers['translation'][0]
            is_pseudo = 'pseudo' in f.qualifiers.keys()

            start_codon_idx = -1 if 'codon_start' not in f.qualifiers.keys() else int(f.qualifiers['codon_start'][0])

            assert((f.type != 'CDS' or translation != '') or is_pseudo)
            assert((f.type == 'CDS' and start_codon_idx != -1) or (f.type != 'CDS'))

            description = f.qualifiers['product']
            gp = GeneProduct(type=f.type,
                             translation=translation,
                             is_pseudo=is_pseudo,
                             start_codon_idx=start_codon_idx,
                             description=description,
                             qualifiers=f.qualifiers)
            all_genes[gene_key].gene_product = gp
    return all_genes


def init_species(record_gb):
    all_genes = init_all_genes(record_gb)
    print(len(all_genes.keys()))

    species_obj = Species(name=record_gb.name,
                          description=record_gb.description,
                          all_genes=all_genes,
                          sequence=str(record_gb.seq))
    return species_obj


def get_stats(record_gb):
    type_stats = {}
    features = record_gb.features
    for i in range(1, len(features)):
        if features[i].type in type_stats:
            type_stats[features[i].type] += 1
        else:
            type_stats[features[i].type] = 1

    return type_stats
