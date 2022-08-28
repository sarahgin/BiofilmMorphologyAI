from enum import Enum


class LEHNINGER_AA_GROUP(Enum):
    POSITIVE = 'P'
    NEGATIVE = 'N'
    AROMATIC = 'R'
    NON_POLAR = 'O'
    POLAR = 'L'


class HYDROPHOBICITY_AA_GROUP(Enum):
    VERY_HYDROPHOBIC = 'V'
    HYDROPHOBIC = 'H'
    NEUTRAL = 'U'
    HYDROPHILIC = 'C'


lehninger_aa_group_list = ['P', 'N', 'R', 'O', 'L']  # TODO


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


lehninger_group_to_aa_dict = {
    LEHNINGER_AA_GROUP.POSITIVE.value: [AA.LYSINE.value, AA.ARGININE.value, AA.HISTIDINE.value],
    LEHNINGER_AA_GROUP.NEGATIVE.value: [AA.GLUTAMATE.value, AA.ASPARTATE.value],
    LEHNINGER_AA_GROUP.AROMATIC.value: [AA.PHENYLALANINE.value, AA.TYROSINE.value, AA.TRYPTOPHANE.value],
    LEHNINGER_AA_GROUP.NON_POLAR.value: [AA.GLYCINE.value, AA.ALANINE.value, AA.PROLINE.value,
                                         AA.VALINE.value, AA.LEUCINE.value,
                                         AA.ISOLEUCINE.value, AA.METHIONINE.value],
    LEHNINGER_AA_GROUP.POLAR.value: [AA.SERINE.value, AA.THREONINE.value, AA.CYSTEINE.value,
                                     AA.ASPARGINE.value, AA.GLUTAMINE.value],
}

# https://www.sigmaaldrich.com/IL/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart
hydrophobicity_group_to_aa_dict = {
    HYDROPHOBICITY_AA_GROUP.VERY_HYDROPHOBIC.value: [AA.ISOLEUCINE.value, AA.LEUCINE.value,
                                                     AA.PHENYLALANINE.value, AA.TRYPTOPHANE.value,
                                                     AA.VALINE.value,  AA.METHIONINE.value],
    HYDROPHOBICITY_AA_GROUP.HYDROPHOBIC.value: [AA.ALANINE.value, AA.TYROSINE.value, AA.CYSTEINE.value],
    HYDROPHOBICITY_AA_GROUP.NEUTRAL.value: [AA.THREONINE.value, AA.GLUTAMATE.value, AA.GLYCINE.value, AA.SERINE.value,
                                            AA.GLUTAMINE, AA.ASPARTATE],
    HYDROPHOBICITY_AA_GROUP.HYDROPHILIC.value: [AA.LYSINE.value, AA.ARGININE.value, AA.HISTIDINE.value, AA.PROLINE.value, AA.ASPARGINE.value],
}
