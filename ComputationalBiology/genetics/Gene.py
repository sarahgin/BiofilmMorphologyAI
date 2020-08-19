from BiofilmMorphologyAI.ComputationalBiology.data_analysis.gene_features_calculator import calculate_feature


class Gene:
    """
    Genomic Sequence on the DNA.
    Attributes:
        DNA sequence (substring of the genome): String (or Seq of biopython)
        Start position
        End position
        Strand # TODO: DNA sequence is different for the minus strand - maybe reverse complement
        Gene Name
        Gene ID
        GeneProduct: represents the gene product (if regulatory: this is None)
    """
    def __init__(self, start: int, end: int, strand, name: str, locus_tag: str): #, sequence, gene_product):
        self.start = start
        self.end = end
        self.strand = strand
        self.locus_tag = locus_tag  # assuming this is unique
        self.name = name
        self.type = type  # TODO: type is enum: one of: 'CDS','rRNA','tRNA','ncRNA','tmRNA', 'regulatory'
        self.id = '{}_{}_{}'.format(start, end, strand)
        #
        # self.gene_product = gene_product  # if type is
        # self.sequence = sequence  #
        # self.gene_features_dict = {}
        # self.protein_features_dict = {}

    def get_id(self) -> str:
        return self.id

    def set_feature(self, feature_name):
        self.gene_features_dict[feature_name] = calculate_feature(feature_name, self.sequence)

    def set_all_features(self, feature_names):
        for f in feature_names:
            self.set_feature(f)

    def is_coding(self):
        return self.gene_product is None
