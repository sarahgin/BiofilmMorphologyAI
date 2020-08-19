from BiofilmMorphologyAI.ComputationalBiology.biologyutils.protein_utils import get_hydrophobic_amino_acids


def compute_hydropobic_aa(seq):
    """
    Given a protein sequence, compute content of hydropobic amino acids
    Return number of matches and percentage of the protein
    input: protein sequence
    output: number of hydropobic amino acids, percentage of hydrophobic amino acids
    """
    aa = get_hydrophobic_amino_acids()
    # TODO
    return 0