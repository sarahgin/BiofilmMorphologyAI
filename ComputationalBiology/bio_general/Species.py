from ComputationalBiology.bio_general.WholeGenomeSequence import WholeGenomeSequence
from ComputationalBiology.bio_general import Gene


class Species:

    def __init__(self, name: str,  description: str, all_genes: dict, sequence: str):
        self.name = name
        self.description = description
        self.all_genes = all_genes
        self.sequence = sequence
        self.genome_sequence = WholeGenomeSequence(sequence)

    def get_gene(self, gene_id: str) -> Gene:
        # TODO: handle if not exists
        # TODO: check if gene_id is unique! in human I saw several gene names of the dame gene
        return self.all_genes[gene_id]

    # def get_gene_sequence(self, gene_id: str) -> str:
    #     """
    #     Given gene-id, get its sequence
    #     """
    #     gene_data = self.all_genes[gene_id]
    #     return self.genome_sequence.get_coding_sequence(gene_data.start, gene_data.end)