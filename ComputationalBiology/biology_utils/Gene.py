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
    def __init__(self, start: int, end: int, strand, type: str, qualifiers, gp=None):

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

    def is_coding(self):
        return self.gene_product is None
