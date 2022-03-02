
class GeneProduct:
    def __init__(self, type, translation, is_pseudo, codon_start, description: str, qualifiers):
        self.type = type
        self.translation = translation
        self.is_pseudo = is_pseudo
        self.codon_start = codon_start
        self.description = description
        self.qualifiers = qualifiers
