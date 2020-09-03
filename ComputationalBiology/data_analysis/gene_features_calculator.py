import re

from ComputationalBiology.bio_general.Gene import Gene
from ComputationalBiology.bio_general.bio_macros import ALL_NT, ValidAlphabet
from ComputationalBiology.bio_general.bio_utils import kmers_generator


def compute_gc_content(dna_sequence: str):
    """
    Given a genomic sequence (str), compute GC content.
    input: DNA sequence (str)
    output:%GC
    """
    return compute_gc_count(dna_sequence) / len(dna_sequence) * 100


def compute_gc_count(dna_sequence: str):
    matches = re.findall(r'[GC]', dna_sequence.upper())
    return len(matches)


def compute_gene_length(dna_sequence: str):
    return len(dna_sequence)


def is_valid_dna(dna_sequence: str):
    matches = re.findall(r'[^ATGC]', dna_sequence.upper())
    return len(matches) == 0


#if not (gene_seq[0:3] == 'ATG' and gene_seq[-3:] in ['TAG', 'TGA', 'TAA']):
#    print((len(gene_seq) - 3), len(protein_seq) * 3)
#    print(lst_CDS[i])


def compute_all_features():
    return 0


def compute_all_kmer_counts(sequence: str, k: int,
                            alphabet: ValidAlphabet,
                            sliding_window=1, start_position=0):
    """

    :param alphabet:
    :param sequence:
    :param k:
    :param sliding_window:
    :param start_position: transcription start site (tss 0 is the +1 site)
    :return:
    """
    # initialize the result:
    result = {}
    gen = kmers_generator(k, alphabet)
    for x in gen:
        result[x] = 0

    # update the result:
    for pos in range(start_position, len(sequence) - k + 1, sliding_window):
        current_kmer = sequence[pos: pos + k]
        result[current_kmer] += 1

    return result


def compute_hexamer_counts(gene: Gene):
    # codon_start is starting from 1
    codon_start = gene.gene_product.codon_start if gene.gene_product is not None else 1

    return compute_all_kmer_counts(gene.coding_sequence,
                                   k=6,
                                   sliding_window=1,
                                   start_position=codon_start-1,  # convert to python indexing
                                   alphabet=ValidAlphabet.NT)


def compute_codon_counts(gene: Gene):
    codon_start = gene.gene_product.codon_start if gene.gene_product is not None else 1

    return compute_all_kmer_counts(gene.coding_sequence,
                                   k=3,
                                   sliding_window=3,
                                   start_position=codon_start-1,
                                   alphabet=ValidAlphabet.NT)


if __name__ == '__main__':

    res = compute_all_kmer_counts('AAACGT',
                            k=3,
                            sliding_window=3,
                            start_position=0,
                            alphabet=ValidAlphabet.NT)
    print(res)

