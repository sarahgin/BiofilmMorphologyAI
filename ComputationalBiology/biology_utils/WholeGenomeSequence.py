# TODO: make this a class with a genome and function to fetch subsequence from the genome sequence
class WholeGenomeSequence:

    def __init__(self, genome_sequence: str):
        self.genome_sequence = genome_sequence
        # self.species_name = species_name
        # self.species_id = species_id

    def get_sub_sequence(self, start_seq: int, end_seq: int):
        """
        " Return the subsequence as given
        """
        return self.genome_sequence[start_seq:end_seq]

    def get_coding_sequence(self, start_seq: int, end_seq: int, strand: int):
        """
        " Return the subsequence as given
        """
        seq = self.get_sub_sequence(start_seq, end_seq)
        if strand == -1:
            seq = seq.reverse_complement()
        return seq

