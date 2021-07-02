import re
import pandas as pd

from ComputationalBiology.bio_general.bio_macros import chem_df

'''
Alanine 	Ala 	A
Arginine 	Arg 	R
Asparagine 	Asn 	N
Aspartic acid 	Asp 	D
Cysteine 	Cys 	C
Glutamic acid 	Glu 	E
Glutamine 	Gln 	Q
Glycine 	Gly 	G
Histidine 	His 	H
Isoleucine 	Ile 	I
Leucine 	Leu 	L
Lysine 	Lys 	K
Methionine 	Met 	M
Phenylalanine 	Phe 	F
Proline 	Pro 	P
Serine 	Ser 	S
Threonine 	Thr 	T
Tryptophan 	Trp 	W
Tyrosine 	Tyr 	Y
Valine 	Val 	V
'''


def compute_hydrophobic_aa(amino_acid_sequence: str) -> float:
    """
    Given a protein sequence, compute content of hydropobic amino acids
    Return number of matches and percentage of the protein
    input: protein sequence
    output: number of hydropobic amino acids, percentage of hydrophobic amino acids
    """
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[GAVLIPFMW]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


def compute_hydrophilic_aa(amino_acid_sequence: str) -> float:
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[KRHDESTCPNQ]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


# AMINO ACID 5 GROUPING
def compute_nonpolar_aa(amino_acid_sequence: str) -> float:
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[GAVLMI]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


def compute_aromatic_aa(amino_acid_sequence: str) -> float:
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[FYW]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


def compute_positive_aa(amino_acid_sequence: str) -> float:
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[KRH]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


def compute_negative_aa(amino_acid_sequence: str) -> float:
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[DE]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


def compute_polar_aa(amino_acid_sequence: str) -> float:
    if len(amino_acid_sequence) == 0:
        return None

    matches = re.findall(r'[STCPNQ]', amino_acid_sequence.upper())
    return len(matches) / len(amino_acid_sequence) * 100


def is_valid_protein(amino_acid_sequence: str) -> bool:
    matches = re.findall(r'[^FSHNGWQTRVLYMCIDAEK]', amino_acid_sequence.upper())
    return len(matches) == 0


def compute_protein_length(amino_acid_sequence: str) -> int:
    return len(amino_acid_sequence)


# chemical properties
def compute_H1(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'H1')


def compute_H2(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'H2')


def compute_H3(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'H3')


def compute_V(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'V')


def compute_P1(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'P1')


def compute_P2(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'P2')


def compute_SASA(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'SASA')


def compute_NCI(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'NCI')


def compute_MASS(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'MASS')


def compute_PKA_COOH(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'PKA_COOH')


def compute_PKA_NH(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'PKA_NH')


def compute_PI(amino_acid_sequence: str):
    return compute_avg_protein_chemical_feature(amino_acid_sequence, 'PI')


def compute_avg_protein_chemical_feature(amino_acid_sequence: str, chemical_feature_name: str):
    s = 0
    for aa in amino_acid_sequence:
        s += float(chem_df[chemical_feature_name][aa])
    return s / len(amino_acid_sequence) if len(amino_acid_sequence) > 0 else 0
