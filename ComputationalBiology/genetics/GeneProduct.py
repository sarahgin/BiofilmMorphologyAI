from BiofilmMorphologyAI.ComputationalBiology.genetics.Gene import Gene


class GeneProduct(Gene):

    # TODO: type is enum: one of: 'CDS','rRNA','tRNA','ncRNA','tmRNA'
    def __init__(self, start, end, sequence, name, id, type, product_description:str):
        self.start = start   # TODO: do we need this?
        self.end = end       # TODO: do we need this?
        self.sequence = sequence
        self.name = name
        self.id = id
        self.type = type
        self.product_description = product_description # str
