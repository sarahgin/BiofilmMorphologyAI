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


aa_group_list = ['P', 'N', 'R', 'O', 'L']


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
                           AA.ASPARGINE.value, AA.GLUTAMINE.value]

}

alphabets_dict = {
    ValidAlphabet.NT: ALL_NT,
    ValidAlphabet.AA: ALL_AA
}

class PTS_EXACT(Enum):
    ER_IMPORT_EXACT = 1
    ER_RETENTION_EXACT = 2
    MITOCHONDRIA_IMPORT_EXACT = 3
    NLS_EXACT = 4
    PEROXISOME_IMPORT_EXACT = 5


protein_targeting_signals_exact_dict = {
    PTS_EXACT.ER_IMPORT_EXACT: '\M\M\S\F\V\S\L\L\L\V\G\I\L\F\W\A\T\E\A\E\Q\L\T\K\C\E\V\F\Q',
    PTS_EXACT.ER_RETENTION_EXACT: '\K\D\E\L',
    PTS_EXACT.MITOCHONDRIA_IMPORT_EXACT: '\M\L\S\L\R\Q\S\I\R\F\F\K\P\A\T\R\T\L\C\S\S\R\Y\L\L',
    PTS_EXACT.NLS_EXACT: '\P\P\K\K\K\R\K\V',
    PTS_EXACT.PEROXISOME_IMPORT_EXACT: '\S\K\L'
}

class PTS_GENERAL(Enum):
    ER_IMPORT = 1
    ER_RETENTION = 2
    NLS = 3
    MITOCHONDRIA_IMPORT = 4
    PEROXISOME_IMPORT = 5


protein_targeting_signals_dict = {
    PTS_GENERAL.ER_IMPORT: '\M\M\S\F\V\SO{7}A{2}O\TN\AN\Q\L\TP\CN\V\F\Q',
    PTS_GENERAL.ER_RETENTION: 'PNN\L',
    PTS_GENERAL.NLS: '\P\A\A\K\R\V\K\L\D',#'\KP.{1}P', #'\P\PP{5}\V',
    PTS_GENERAL.MITOCHONDRIA_IMPORT: '\M\L\S\LP\Q\S\IP\F\FP\P\A\TP\T\L\C\S\SP\Y\L\L',
    PTS_GENERAL.PEROXISOME_IMPORT: '\SP\L'
}
