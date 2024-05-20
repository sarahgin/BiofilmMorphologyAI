import re

from BioinformaticsLab.ComputationalBiology.bio_general.Gene import Gene


def create_prefix_suffix_dict(seq: str,
                              prefix_length_min,
                              prefix_length_max,
                              suffix_length_min,
                              suffix_length_max,
                              codon_start=0):
    full_dict = {}
    counts_dict = {}
    relative_positions_dict = {}
    seq_length = len(seq)

    # Note: if there are UTRs before 'AUG' (codon start), they are currently ignored
    # gene.codon_start might be ~= 0
    for pos in range(codon_start, len(seq) - 1):
        for curr_prefix_len in range(prefix_length_min, prefix_length_max + 1):
            curr_prefix = seq[pos: pos + curr_prefix_len]
            for curr_suffix_len in range(suffix_length_min, suffix_length_max + 1):
                if pos + curr_prefix_len + curr_suffix_len > len(seq):
                    continue
                curr_suffix = seq[pos + curr_prefix_len: pos + curr_prefix_len + curr_suffix_len]

                if curr_prefix not in full_dict.keys():
                    full_dict[curr_prefix] = {}
                    #assert (curr_prefix not in prefix_counts_dict.keys())
                    counts_dict[curr_prefix] = 1
                    relative_positions_dict[curr_prefix] = [pos / seq_length]
                else:
                    counts_dict[curr_prefix] += 1
                    relative_positions_dict[curr_prefix].append(pos / seq_length)

                if curr_suffix not in full_dict[curr_prefix].keys():
                    full_dict[curr_prefix][curr_suffix] = 1
                else:
                    full_dict[curr_prefix][curr_suffix] += 1

    return full_dict, counts_dict, relative_positions_dict


def create_gene_prefix_suffix_dict(gene: Gene,
                                   prefix_length_min=None,
                                   prefix_length_max=None,
                                   suffix_length_min=None,
                                   suffix_length_max=None):
    gene_length = len(gene.coding_sequence)
    if prefix_length_min is None:
        prefix_length_min = 4
    if prefix_length_max is None:
        prefix_length_max = 4  # gene_length
    if suffix_length_min is None:
        suffix_length_min = 2
    if suffix_length_max is None:
        suffix_length_max = 2  # gene_length

    return create_prefix_suffix_dict(gene.coding_sequence,
                                          prefix_length_min,
                                          prefix_length_max,
                                          suffix_length_min,
                                          suffix_length_max,
                                          gene.codon_start)


def create_protein_prefix_suffix_dict(gene: Gene,
                                   prefix_length_min,
                                   prefix_length_max,
                                   suffix_length_min,
                                   suffix_length_max):

    protein_length = gene.gene_product.translation
    if prefix_length_min is None:
        prefix_length_min = 1
    if prefix_length_max is None:
        prefix_length_max = protein_length
    if suffix_length_min is None:
        suffix_length_min = 1
    if suffix_length_max is None:
        suffix_length_max = protein_length
    return create_prefix_suffix_dict(gene.gene_product.translation,
                                          prefix_length_min,
                                          prefix_length_max,
                                          suffix_length_min,
                                          suffix_length_max,
                                          gene.codon_start)


if __name__ == '__main__':
    full_dict, counts_dict, relative_positions_dict = create_prefix_suffix_dict(seq='AGTCTAGAGTCTADDDDDDD',
                                                            prefix_length_min=2,
                                                            prefix_length_max=2,
                                                            suffix_length_min=3,
                                                            suffix_length_max=3,
                                                            codon_start=0)

    print(full_dict)
    print(counts_dict)
    print(relative_positions_dict)
    print('Done')