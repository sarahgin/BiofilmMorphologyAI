from ComputationalBiology.biology_utils.WholeGenomeSequence import WholeGenomeSequence
from ComputationalBiology.biology_utils import Gene


class Species:

    def __init__(self, species_name: str,  species_id: str, all_genes: dict, genome_seq: str):
        self.species_name = species_name
        self.species_id = species_id
        self.genome_sequence = WholeGenomeSequence(genome_seq)

        # Collection of genes: dictionary of string-unique id(key)-->Gene(value)
        self.all_genes = all_genes

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