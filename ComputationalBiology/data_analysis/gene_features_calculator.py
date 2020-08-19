#  Compute different feature for a Gene object

def calculate_feature(feature_name, seq):
    if feature_name == 'gc_content':
        return compute_gc_content(seq)
    else:
        return None #  ivnalid feature name

def compute_gc_content(seq):
    """
    Given a genomic sequence (String or Seq), compute GC content.
    Return #of G+C, #of total nucleotides
    input: DNA sequence (str or Seq)
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
    # TODO: add  warning if needed, e.g. illegal nucleotide
    return counter_C + counter_G,  (counter_C + counter_G) / len(seq)

