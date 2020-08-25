from ComputationalBiology.biology_utils.Gene import Gene


class GeneProduct:
    def __init__(self, type, translation, is_pseudo, start_codon_idx, qualifiers):
        self.type = type
        self.translation = translation
        self.is_pseudo = is_pseudo
        self.start_codon_idx = start_codon_idx
        self.qualifiers = qualifiers

