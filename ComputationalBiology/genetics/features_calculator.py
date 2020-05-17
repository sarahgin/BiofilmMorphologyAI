from BiofilmMorphologyAI.ComputationalBiology.biologyutils.GenomeSequence import GenomeSequence
from BiofilmMorphologyAI.ComputationalBiology.biologyutils.protein_utils import get_hydrophobic_amino_acids

# TODO: make a class of features calculator?
# class FeaturesCalculator:
#     def __init__(self, genome_seq):
#         self.genome = Genome(genome_seq)

#
# def get_single_feature(sequence_obj, feature):
#     pass

def get_tuple_features(single_tuple, genome):
    """
    Returns a list of features, currently supports:
    1. GC content
    2. % of hydropobic amino-acids
    """
    # fetch the gene info:# for 1st gene:
    start_s = single_tuple.gene.location.start.position
    end_s = single_tuple.gene.location.end.position
    strand = single_tuple.gene.location.strand

    gene_sequence = genome.get_reference_seq(start_s, end_s, strand)  # string
    print(gene_sequence)
    # prot_sequence = t.gene_product

    gene_features = get_all_gene_features(gene_sequence)
    # prot_features = get_prot_features(prot_sequence)
    all_features = {}
    all_features.update(gene_features)
    # all_features.update(prot_features)
    return all_features


def get_all_gene_features(gene_seq):
    """
    Computes all features related to the gene sequence (DNA)
    """
    gc_count, gc_percent = compute_gc_content(gene_seq)
    result = {'gc_count': gc_count,  'gc_content': gc_percent}
    return result


def get_protein_features(prot_seq):
    """
    Computes all features related to the protein sequence
    """
    hydropobic_aa_count, hydropobic_aa_percent = compute_hydropobic_aa(prot_seq)
    result = {'hydropobic_aa_percent': hydropobic_aa_percent}
    return result



def compute_gc_content(seq):
    """
    Given a genomic sequence, compute GC content. Return #of G+C, #of total nucleotides
    input: DNA sequence
    output:#GC, #GC/(#A+#T+#C+#G)
    """
    gene_len = len(seq)
    counter_G = seq.count('G')
    counter_C = seq.count('C')
    counter_A = seq.count('A')
    counter_T = seq.count('T')
    print('counter_A:{}, counter_C:{}, counter_G:{}, counter_T: {}'.format(
        counter_A, counter_C, counter_G, counter_T))
    sum_nt = sum([counter_A, counter_C, counter_G, counter_T])
    print('sum_nt:', sum_nt)
    print('gene_len:', gene_len)
    # assert (gene_len == sum_nt), 'genome_len: {}, sum_nt: {} '.format(gene_len, sum_nt)
    return counter_C + counter_G,  (counter_C + counter_G) / len(seq)


def compute_hydropobic_aa(prot_seq):
    """
    Given a protein sequence, compute content of hydropobic amino acids
    Return number of matches and percentage of the protein
    input: protein sequence
    output: number of hydropobic amino acids,  percentage of hydropobic amino acids
    """
    aa = get_hydrophobic_amino_acids()
    pass



