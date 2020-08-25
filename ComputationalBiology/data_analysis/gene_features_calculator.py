import re


def compute_gc_content(dna_sequence):
    """
    Given a genomic sequence (str), compute GC content.
    input: DNA sequence (str)
    output:%GC
    """
    matches = re.findall(r'[GC]', dna_sequence.upper())
    return len(matches), len(matches) / len(dna_sequence) * 100


def compute_melting_point(dna_sequence):
    #Tm = 4(G + C) + 2(A + T)
    gc_count, _ = compute_gc_content(dna_sequence)
    return 4*gc_count + 2*(compute_length(dna_sequence) - gc_count)


def compute_length(dna_sequence):
    return len(dna_sequence)


def is_valid_dna(dna_sequence):
    matches = re.findall(r'[^ATGC]', dna_sequence.upper())
    return len(matches) == 0


#if not (gene_seq[0:3] == 'ATG' and gene_seq[-3:] in ['TAG', 'TGA', 'TAA']):
#    print((len(gene_seq) - 3), len(protein_seq) * 3)
#    print(lst_CDS[i])


def compute_all_features():
    return 0


if __name__ == '__main__':
    print(is_valid_dna('ATAG'))

