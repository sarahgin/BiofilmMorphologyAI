import re
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
    matches = re.findall(r'[GAVLIPFMW]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def compute_hydrophilic_aa(amino_acid_sequence: str) -> float:
    matches = re.findall(r'[KRHDESTCPNQ]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


#AMINO ACID 5 GROUPING
def compute_nonpolar_aa(amino_acid_sequence: str) -> float:
    matches = re.findall(r'[GAVLMI]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def compute_aromatic_aa(amino_acid_sequence: str) -> float:
    matches = re.findall(r'[FYW]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def compute_positive_aa(amino_acid_sequence: str) -> float:
    matches = re.findall(r'[KRH]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def compute_negative_aa(amino_acid_sequence: str) -> float:
    matches = re.findall(r'[DE]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def compute_polar_aa(amino_acid_sequence: str) -> float:
    matches = re.findall(r'[STCPNQ]', amino_acid_sequence.upper())
    return len(matches)/len(amino_acid_sequence)*100


def is_valid_protein(amino_acid_sequence: str) -> bool:
    matches = re.findall(r'[^FSHNGWQTRVLYMCIDAEK]', amino_acid_sequence.upper())
    return len(matches) == 0


def compute_protein_length(amino_acid_sequence: str) -> int:
    return len(amino_acid_sequence)


def compute_all_features():
    return 0

if __name__ == '__main__':
    print(is_valid_protein('AQAX'))
