import re


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


if __name__ == '__main__':
    print(is_valid_dna('ATAG'))

