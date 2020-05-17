# TODO: make this a class with a genome and function to fetch subsequence from the genome sequence
class GenomeSequence:
    # sepecies
    # id = 3610

    def __init__(self):
        self.genome_sequence = ''

    def __init__(self, genome_sequence):
        self.genome_sequence = genome_sequence
        # self.id = id

    def get_reference_seq(self, start_seq, end_seq):
        return self.genome_sequence[start_seq:end_seq]

    # return the sequence
    def get_coding_sequence(self, start_seq, end_seq, strand):
        seq = self.get_reference_seq(start_seq, end_seq)
        if strand == -1:
            seq = seq.reverse_complement()
        return seq

