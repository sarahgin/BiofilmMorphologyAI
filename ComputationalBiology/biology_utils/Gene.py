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
    def __init__(self, start: int, end: int, strand, type: str, coding_sequence, qualifiers, gp=None):

        self.start = start
        self.end = end
        self.strand = strand
        self.type = type
        self.coding_sequence = coding_sequence
        self.id = self.create_id()

        #raw data
        self.qualifiers = qualifiers
        self.gene_product = gp

    def create_id(self):
        return '{}_{}_{}'.format(self.start, self.end, self.strand)

    def get_id(self) -> str:
        return self.id

    def has_product(self):
        return self.gene_product is None
