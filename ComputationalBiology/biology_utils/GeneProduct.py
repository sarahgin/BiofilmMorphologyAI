from ComputationalBiology.biology_utils.Gene import Gene


class GeneProduct(Gene):
    def __init__(self, type, translation, is_pseudo, qualifiers):
        self.type = type
        self.translation = translation
        self.is_pseudo = is_pseudo
        self.qualifiers = qualifiers

