import re


def compute_hydrophobic_aa(amino_acid_sequence):
    """
    Given a protein sequence, compute content of hydropobic amino acids
    Return number of matches and percentage of the protein
    input: protein sequence
    output: number of hydropobic amino acids, percentage of hydrophobic amino acids
    """
    matches = re.findall(r'[GAVLIPFMW]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def is_valid_protein(amino_acid_sequence):
    matches = re.findall(r'[^FSHNGWQTRVLYMCIDAEK]', amino_acid_sequence.upper())
    return len(matches) == 0


def compute_length(amino_acid_sequence):
    return len(amino_acid_sequence)


def compute_all_features():
    return 0

if __name__ == '__main__':
    print(is_valid_protein('AQAX'))
