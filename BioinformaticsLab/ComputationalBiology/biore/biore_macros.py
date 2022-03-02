from enum import Enum

genetic_code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

ALL_NT = 'ACGT'

ALL_AA = 'GAVLIPFMWKRHDESTCPNQ'


class ValidAlphabet(Enum):
    NT = 1
    AA = 2


class AA_GROUP(Enum):
    POSITIVE = 'P'
    NEGATIVE = 'N'
    AROMATIC = 'R'
    NON_POLAR = 'O'
    POLAR = 'L'
    HYDROPHOBIC = 'B'
    HYDROPHILIC = 'Y'


class HYDROPHOBICITY_AA_GROUP(Enum):
    VERY_HYDROPHOBIC = 'P'
    HYDROPHOBIC = 'N'
    NEUTRAL = 'R'
    HYDROPHILIC = 'Q'


aa_group_list = ['P', 'N', 'R', 'O', 'L', 'B', 'Y']

class NT(Enum):
    ADENINE = 'A'
    CYTOSINE = 'C'
    GUANINE = 'G'
    THYMINE = 'T'

class AA(Enum):
    ALANINE = 'A'
    CYSTEINE = 'C'
    ASPARTATE = 'D'
    GLUTAMATE = 'E'
    PHENYLALANINE = 'F'
    GLYCINE = 'G'
    HISTIDINE = 'H'
    ISOLEUCINE = 'I'
    LYSINE = 'K'
    LEUCINE = 'L'
    METHIONINE = 'M'
    ASPARGINE = 'N'
    PROLINE = 'P'
    GLUTAMINE = 'Q'
    ARGININE = 'R'
    SERINE = 'S'
    THREONINE = 'T'
    VALINE = 'V'
    TRYPTOPHANE = 'W'
    TYROSINE = 'Y'


aa_to_BY_dict = {AA.ALANINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.CYSTEINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.ASPARTATE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.GLUTAMATE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.PHENYLALANINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.GLYCINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.HISTIDINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.ISOLEUCINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.LYSINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.LEUCINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.METHIONINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.ASPARGINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.PROLINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.GLUTAMINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.ARGININE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.SERINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.THREONINE.value: AA_GROUP.HYDROPHILIC.value,
                    AA.VALINE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.TRYPTOPHANE.value: AA_GROUP.HYDROPHOBIC.value,
                    AA.TYROSINE.value: AA_GROUP.HYDROPHOBIC.value
                    }

aa_to_group_dict = {AA.ALANINE.value: AA_GROUP.NON_POLAR.value,
                    AA.CYSTEINE.value: AA_GROUP.POLAR.value,
                    AA.ASPARTATE.value: AA_GROUP.NEGATIVE.value,
                    AA.GLUTAMATE.value: AA_GROUP.NEGATIVE.value,
                    AA.PHENYLALANINE.value: AA_GROUP.AROMATIC.value,
                    AA.GLYCINE.value: AA_GROUP.NON_POLAR.value,
                    AA.HISTIDINE.value: AA_GROUP.POSITIVE.value,
                    AA.ISOLEUCINE.value: AA_GROUP.NON_POLAR.value,
                    AA.LYSINE.value: AA_GROUP.POSITIVE.value,
                    AA.LEUCINE.value: AA_GROUP.NON_POLAR.value,
                    AA.METHIONINE.value: AA_GROUP.NON_POLAR.value,
                    AA.ASPARGINE.value: AA_GROUP.POLAR.value,
                    AA.PROLINE.value: AA_GROUP.NON_POLAR.value,
                    AA.GLUTAMINE.value: AA_GROUP.POLAR.value,
                    AA.ARGININE.value: AA_GROUP.POSITIVE.value,
                    AA.SERINE.value: AA_GROUP.POLAR.value,
                    AA.THREONINE.value: AA_GROUP.POLAR.value,
                    AA.VALINE.value: AA_GROUP.NON_POLAR.value,
                    AA.TRYPTOPHANE.value: AA_GROUP.AROMATIC.value,
                    AA.TYROSINE.value: AA_GROUP.AROMATIC.value
                    }

aa_to_hydrophobicity_group_dict = {
                    AA.ALANINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHOBIC.value,
                    AA.CYSTEINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHOBIC.value,
                    AA.ASPARTATE.value: HYDROPHOBICITY_AA_GROUP.NEUTRAL.value,
                    AA.GLUTAMATE.value: HYDROPHOBICITY_AA_GROUP.NEUTRAL.value,
                    AA.PHENYLALANINE.value: HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value,
                    AA.GLYCINE.value: HYDROPHOBICITY_AA_GROUP.NEUTRAL.value,
                    AA.HISTIDINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value,
                    AA.ISOLEUCINE.value: HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value,
                    AA.LYSINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value,
                    AA.LEUCINE.value: HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value,
                    AA.METHIONINE.value: HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value,
                    AA.ASPARGINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value,
                    AA.PROLINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value,
                    AA.GLUTAMINE.value: HYDROPHOBICITY_AA_GROUP.NEUTRAL.value,
                    AA.ARGININE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value,
                    AA.SERINE.value: HYDROPHOBICITY_AA_GROUP.NEUTRAL.value,
                    AA.THREONINE.value: HYDROPHOBICITY_AA_GROUP.NEUTRAL.value,
                    AA.VALINE.value: HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value,
                    AA.TRYPTOPHANE.value: HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value,
                    AA.TYROSINE.value: HYDROPHOBICITY_AA_GROUP.HYDROPHOBIC.value
                    }

group_to_aa_dict = {
    AA_GROUP.POSITIVE.value: [AA.LYSINE.value, AA.ARGININE.value, AA.HISTIDINE.value],
    AA_GROUP.NEGATIVE.value: [AA.GLUTAMATE.value, AA.ASPARTATE.value],
    AA_GROUP.AROMATIC.value: [AA.PHENYLALANINE.value, AA.TYROSINE.value, AA.TRYPTOPHANE.value],
    AA_GROUP.NON_POLAR.value: [AA.GLYCINE.value, AA.ALANINE.value, AA.PROLINE.value,
                               AA.VALINE.value, AA.LEUCINE.value,
                               AA.ISOLEUCINE.value, AA.METHIONINE.value],
    AA_GROUP.POLAR.value: [AA.SERINE.value, AA.THREONINE.value, AA.CYSTEINE.value,
                           AA.ASPARGINE.value, AA.GLUTAMINE.value],
}

# https://www.sigmaaldrich.com/IL/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart
HYDROPHOBICITY_group_to_aa_dict = {
    HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value: [AA.ISOLEUCINE.value, AA.LEUCINE.value,
                                                      AA.PHENYLALANINE.value, AA.TRYPTOPHANE.value,
                                                      AA.VALINE.value,  AA.METHIONINE.value],
    HYDROPHOBICITY_AA_GROUP.HYDROPHOBIC.value: [AA.ALANINE.value,AA.TYROSINE.value, AA.CYSTEINE.value],
    HYDROPHOBICITY_AA_GROUP.NEUTRAL.value: [AA.THREONINE.value, AA.GLUTAMATE.value, AA.GLYCINE.value, AA.SERINE.value,
                                            AA.GLUTAMINE, AA.ASPARTATE],
    HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value: [AA.LYSINE.value, AA.ARGININE.value, AA.HISTIDINE.value, AA.PROLINE.value, AA.ASPARGINE.value],
}

alphabets_dict = {
    ValidAlphabet.NT: ALL_NT,
    ValidAlphabet.AA: ALL_AA
}


class PTS(Enum):
    ER_IMPORT_1 = '\L\L\L\V\G\I\L\F\W\A.\E.\E.{3}\K.\E'  # Lehninger exact
    ER_IMPORT_2 = 'B{10}.N.N.{3}P.N'
    ER_IMPORT_3 = 'P+B{6-12}'  # pts5
    ER_IMPORT_ALL = '(' + ER_IMPORT_1 + ')|' + \
                    '(' + ER_IMPORT_2 + ')|' + \
                    '(' + ER_IMPORT_3 + ')'

    MITOCHONDRIA_IMPORT_0 = '[^\R\K]{4}\R[^\R\K]{3}\R[^\R\K]{2}\K[^\R\K]{3}\R[^\R\K]{5}\R[^\R\K]{3}'  # Lehninger exact
    MITOCHONDRIA_IMPORT_1 = '[^\R\K]{4}[\R\K][^\R\K]{3}[\R\K][^\R\K]{2}[\R\K][^\R\K]{3}[\R\K][^\R\K]{5}[\R\K][^\R\K]{3}'  # Lehninger exact
    MITOCHONDRIA_IMPORT_2 = '[^P]{4}P[^P]{3}P[^P]{2}P[^P]{3}P[^P]{5}P[^P]{3}'  # Lehninger general
    MITOCHONDRIA_IMPORT_3 = '([^\E\\N\R\K][\R\K]){3,5}'  # pts5
    MITOCHONDRIA_IMPORT_4 = '([^\R\K]{2,5}[\R\K]){3,5}'
    MITOCHONDRIA_IMPORT_ALL = '(' + MITOCHONDRIA_IMPORT_0 + ')|' + \
                              '(' + MITOCHONDRIA_IMPORT_1 + ')|' + \
                              '(' + MITOCHONDRIA_IMPORT_2 + ')|' + \
                              '(' + MITOCHONDRIA_IMPORT_3 + ')|' + \
                              '(' + MITOCHONDRIA_IMPORT_4 + ')'

    NLS_1 = '\K\K\K\R\K'
    NLS_2 = '\P\A\A\K\R\V\K\L\D'
    NLS_3 = '\KP.{1}P'
    NLS_4 = '\P\PP{5}\V'
    NLS_5 = 'P{5}'
    NLS_6 = 'P{2,4}.{9,11}P{2,4}'  # pts5
    NLS_7 = 'P{2}.{10}P{2}'
    NLS_ALL = '(' + NLS_1 + ')|' + \
              '(' + NLS_2 + ')|' + \
              '(' + NLS_3 + ')|' + \
              '(' + NLS_4 + ')|' + \
              '(' + NLS_5 + ')|' + \
              '(' + NLS_6 + ')|' + \
              '(' + NLS_7 + ')'

    # try later maybe
    ER_RETENTION_1 = '\K\D\E\L'  # Lehninger exact
    ER_RETENTION_2 = 'PNN\L'  # Lehninger general
    PEROXISOME_IMPORT_1 = '\S\K\L'
    PEROXISOME_IMPORT_2 = '\SP\L'

    # alpha helix regexes
    ALPHA_HELIX_1 = '[\A\E\L\M\K]{10}'
    ALPHA_HELIX_2 = '[^\P\G\Y\S]{10}'
    ALPHA_HELIX_ALL = '(' + ALPHA_HELIX_1 + ')|' + \
                      '(' + ALPHA_HELIX_2 + ')|'


pts_to_run = [PTS.NLS_ALL]
# pts_to_run = [PTS.NLS_1, PTS.NLS_2, PTS.NLS_3, PTS.NLS_4, PTS.NLS_5, PTS.NLS_6, PTS.NLS_ALL]
# pts_to_run = [PTS.MITOCHONDRIA_IMPORT_0, PTS.MITOCHONDRIA_IMPORT_1, PTS.MITOCHONDRIA_IMPORT_2, PTS.MITOCHONDRIA_IMPORT_3, PTS.MITOCHONDRIA_IMPORT_4, PTS.MITOCHONDRIA_IMPORT_ALL]
# pts_to_run = [PTS.ER_IMPORT_1, PTS.ER_IMPORT_2, PTS.ER_IMPORT_3, PTS.ER_IMPORT_ALL]
# pts_to_run = [PTS.ALPHA_HELIX_1, PTS.ALPHA_HELIX_2, PTS.ALPHA_HELIX_ALL]
