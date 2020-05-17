
class DNASequence:
    """
    Genomic Sequence on the DNA.
    Attributes:
        DNA sequence (substring of the genome)
        Start position
        End position
        Strand # TODO: DNA sequence is different for the minus strand - maybe reverse complement
        Sequence ID
        Type: "gene" or "regulatory"
        DNASequenceProduct: represents the gene product (if regulatory: this in null)
    """
    def __init__(self, genbank_tuple, genome_object):
        gene = genbank_tuple.gene
        self.sequence = genome_object.get_sub_seq()  # get sequence from tuple
        genbank_tuple.gene # start, end, stand
        genbank_tuple.gene_product
        # start_s = features[1].location.start.position
        # end_s = features[1].location.end.position
        # strand = features[1].location.strand
        pass

    def __init__(self, sequence, type):
        self.sequence = sequence
        self.type = type
        pass

    def get_feature(self, feature_name):
        # TODO: feature name is enum of all features
        if feature_name == 'gc_content':
            return get_gc_content(self)