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


def compute_a_content(dna_sequence: str):
    return compute_a_count(dna_sequence) / len(dna_sequence) * 100


def compute_a_count(dna_sequence: str):
    matches = re.findall(r'[A]', dna_sequence.upper())
    return len(matches)


def compute_gc_count(dna_sequence: str):
    matches = re.findall(r'[GC]', dna_sequence.upper())
    return len(matches)


def compute_gene_length(dna_sequence: str):
    return len(dna_sequence)


def is_valid_dna(dna_sequence: str):
    matches = re.findall(r'[^ATGC]', dna_sequence.upper())
    return len(matches) == 0


# if not (gene_seq[0:3] == 'ATG' and gene_seq[-3:] in ['TAG', 'TGA', 'TAA']):
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


def compute_all_kmer_positions(sequence: str, k: int,
                               sliding_window=1, start_position=0):
    result = {}
    for pos in range(start_position, len(sequence) - k + 1, sliding_window):
        current_kmer = sequence[pos: pos + k]
        if current_kmer in result.keys():
            result[current_kmer].append(pos)
        else:
            result[current_kmer] = [pos]

    return result



def create_prefix_suffix_dict(seq: str,
                              prefix_length_min,
                              prefix_length_max,
                              suffix_length_min,
                              suffix_length_max,
                              codon_start=0):
    # seq = gene.coding_sequence
    result = {}
    prefix_counts_dict = {}
    # Note: if there are UTRs before 'AUG' (codon start), they are currently ignored
    # gene.codon_start might be ~= 0
    for pos in range(codon_start, len(seq) - 1):
        for curr_prefix_len in range(prefix_length_min, prefix_length_max + 1):
            curr_prefix = seq[pos: pos + curr_prefix_len]
            for curr_suffix_len in range(suffix_length_min, suffix_length_max + 1):
                if pos + curr_prefix_len + curr_suffix_len > len(seq):
                    continue
                curr_suffix = seq[pos + curr_prefix_len: pos + curr_prefix_len + curr_suffix_len]

                if curr_prefix not in result.keys():
                    result[curr_prefix] = {}
                    assert (curr_prefix not in prefix_counts_dict.keys())
                    prefix_counts_dict[curr_prefix] = 1
                else:
                    prefix_counts_dict[curr_prefix] += 1

                if curr_suffix not in result[curr_prefix].keys():
                    result[curr_prefix][curr_suffix] = 1
                else:
                    result[curr_prefix][curr_suffix] += 1
    return result, prefix_counts_dict


def create_gene_prefix_suffix_dict(gene: Gene,
                                   prefix_length_min,
                                   prefix_length_max,
                                   suffix_length_min,
                                   suffix_length_max):

    return create_prefix_suffix_dict(gene.coding_sequence,
                                          prefix_length_min,
                                          prefix_length_max,
                                          suffix_length_min,
                                          suffix_length_max,
                                          gene.codon_start)


def compute_hexamer_positions(gene: Gene):
    # codon_start is starting from 1
    codon_start = gene.gene_product.codon_start if gene.gene_product is not None else 1

    return compute_all_kmer_positions(gene.coding_sequence,
                                      k=6,
                                      sliding_window=1,
                                      start_position=codon_start - 1)  # convert to python indexing


def compute_hexamer_next_nucleotide(gene: Gene):
    codon_start = gene.gene_product.codon_start if gene.gene_product is not None else 1
    return compute_all_suffixes_given_prefix(gene.coding_sequence,
                                             prefix_length=6,
                                             suffix_length=4,
                                             start_position=codon_start - 1)

#start_position is the first position from which we begin the prefix analysis.
#Use **0** when starting from the first nucleotide of the sequence
#Use **gene.codon_start** when staring the sequence from the first amino-acid coding codon
#def compute_all_suffixes_given_prefix(sequence: str,
#                                      prefix_length: int,
#                                      suffix_length: int,
#                                      start_position: int):
#    result = {}
#    for pos in range(start_position, len(sequence) - prefix_length + 1 - suffix_length):
#        current_kmer = sequence[pos: pos + prefix_length]
#        next_nucleotide = sequence[pos + prefix_length: pos + prefix_length + suffix_length]
#        if current_kmer not in result.keys():
#            result[current_kmer] = {}
#        if next_nucleotide not in result[current_kmer].keys():
#            result[current_kmer][next_nucleotide] = 0
#
#        result[current_kmer][next_nucleotide] += 1
#    return result


def compute_codon_counts(gene: Gene):
    return compute_all_kmer_counts(gene.coding_sequence,
                                   k=3,
                                   sliding_window=3,
                                   start_position=gene.codon_start,
                                   alphabet=ValidAlphabet.NT)


if __name__ == '__main__':
    res = compute_all_kmer_counts('AAACGT',
                                  k=3,
                                  sliding_window=3,
                                  start_position=0,
                                  alphabet=ValidAlphabet.NT)
    print(res)
