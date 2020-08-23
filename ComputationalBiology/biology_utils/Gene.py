from ComputationalBiology.data_analysis.gene_features_calculator import calculate_feature


class Gene:
    """
    Genomic Sequence on the DNA.
    Attributes:
        Start position
        End position
        Strand
        Type - gene or regulatory
        GeneProduct: represents the gene product
    """
    def __init__(self, start: int, end: int, strand, type: str, qualifiers, gp = None):

        self.start = start
        self.end = end
        self.strand = strand
        self.type = type
        self.id = self.create_id()

        #raw data
        self.qualifiers = qualifiers
        self.gene_product = gp

    def create_id(self):
        return '{}_{}_{}'.format(self.start, self.end, self.strand)

    def get_id(self) -> str:
        return self.id

    def set_feature(self, feature_name):
        self.gene_features_dict[feature_name] = calculate_feature(feature_name, self.sequence)

    def set_all_features(self, feature_names):
        for f in feature_names:
            self.set_feature(f)

    def is_coding(self):
        return self.gene_product is None
