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


aa_group_list = ['P', 'N', 'R', 'O', 'L', 'B', 'Y']


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


group_to_amino_acid_dict = {
    AA_GROUP.POSITIVE.value: [AA.LYSINE.value, AA.ARGININE.value, AA.HISTIDINE.value],
    AA_GROUP.NEGATIVE.value: [AA.GLUTAMATE.value, AA.ASPARTATE.value],
    AA_GROUP.AROMATIC.value: [AA.PHENYLALANINE.value, AA.TYROSINE.value, AA.TRYPTOPHANE.value],
    AA_GROUP.NON_POLAR.value: [AA.GLYCINE.value, AA.ALANINE.value, AA.PROLINE.value,
                               AA.VALINE.value, AA.LEUCINE.value,
                               AA.ISOLEUCINE.value, AA.METHIONINE.value],
    AA_GROUP.POLAR.value: [AA.SERINE.value, AA.THREONINE.value, AA.CYSTEINE.value,
                           AA.ASPARGINE.value, AA.GLUTAMINE.value],
    AA_GROUP.HYDROPHOBIC.value: [AA.ALANINE.value, AA.VALINE.value, AA.ISOLEUCINE.value, AA.LEUCINE.value,
                           AA.METHIONINE.value, AA.PHENYLALANINE.value, AA.TYROSINE.value, AA.TRYPTOPHANE.value],
    AA_GROUP.HYDROPHILIC.value: [AA.LYSINE.value, AA.ARGININE.value, AA.HISTIDINE.value,
                           AA.ASPARTATE.value, AA.GLUTAMATE.value,
                           AA.SERINE.value, AA.THREONINE.value, AA.ASPARGINE.value, AA.GLUTAMINE.value]

}

alphabets_dict = {
    ValidAlphabet.NT: ALL_NT,
    ValidAlphabet.AA: ALL_AA
}


class PTS(Enum):
    ER_IMPORT_1 = '\M\M\S\F\V\S\L\L\L\V\G\I\L\F\W\A\T\E\A\E\Q\L\T\K\C\E\V\F\Q'  # Lehninger exact
    ER_IMPORT_2 = '\M\M\S\F\V\SO{7}A{2}O\TN\AN\Q\L\TP\CN\V\F\Q'  # Lehninger general
    ER_IMPORT_3 = 'P+B{6-12}'  # pts5

    ER_RETENTION_1 = '\K\D\E\L'  # Lehninger exact
    ER_RETENTION_2 = 'PNN\L'  # Lehninger general

    MITOCHONDRIA_IMPORT_1 = '\M\L\S\L\R\Q\S\I\R\F\F\K\P\A\T\R\T\L\C\S\S\R\Y\L\L'  # Lehninger exact
    MITOCHONDRIA_IMPORT_2 = '\M\L\S\LP\Q\S\IP\F\FP\P\A\TP\T\L\C\S\SP\Y\L\L'  # Lehninger general
    MITOCHONDRIA_IMPORT_3 = r'([^\E\N\R\K][\R\K]){3,5}'  # pts5

    NLS_1 = '\P\P\K\K\K\R\K\V'
    NLS_2 = '\P\A\A\K\R\V\K\L\D'
    NLS_3 = '\KP.{1}P'
    NLS_4 = '\P\PP{5}\V'
    NLS_5 = '(P{5})|(P{2,4}.{10}P{2,4})'  # pts5

    PEROXISOME_IMPORT_1 = '\S\K\L'
    PEROXISOME_IMPORT_2 = '\SP\L'


pts_to_run = {PTS.ER_IMPORT_3, PTS.MITOCHONDRIA_IMPORT_3, PTS.NLS_5}
