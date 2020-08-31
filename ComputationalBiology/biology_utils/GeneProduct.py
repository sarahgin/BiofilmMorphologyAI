
class GeneProduct:
    def __init__(self, type, translation, is_pseudo, start_codon_idx, description: str, qualifiers):
        self.type = type
        self.translation = translation
        self.is_pseudo = is_pseudo
        self.start_codon_idx = start_codon_idx
        self.description = description
        self.qualifiers = qualifiers
