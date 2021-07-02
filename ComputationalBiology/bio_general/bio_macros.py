from enum import Enum
import pandas as pd
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

alphabets_dict = {
    ValidAlphabet.NT: ALL_NT,
    ValidAlphabet.AA: ALL_AA
}

PREFIX_LENGTH_MIN = 5
PREFIX_LENGTH_MAX = 5
SUFFIX_LENGTH_MIN = 5
SUFFIX_LENGTH_MAX = 5


chemicals_features = ['H1', 'H2', 'H3', 'V', 'P1', 'P2', 'SASA', 'NCI', 'MASS', 'PKA_COOH', 'PKA_NH', 'PI']
chemical_properties = {'A': ['0.62', '-0.5', '2', '27.5', '8.1', '0.046', '1.181', '0.007187', '71.0788', '2.35', '9.87', '6.11'],
	                   'C': ['0.29', '-1', '2', '44.6', '5.5', '0.128', '1.461', '-0.03661', '103.1388', '1.92', '10.70', '5.15'],
	                   'D': ['-0.9', '3', '4', '40', '13', '0.105', '1.587', '-0.02382', '115.0886', '1.99', '9.90', '2.98'],
	                   'E': ['-0.74', '3', '4', '62', '12.3', '0.151', '1.862', '0.006802', '129.1155', '2.10', '9.47', '3.08'],
	                   'F': ['1.19', '-2.5', '2', '115.5', '5.2', '0.29', '2.228', '0.037552', '147.1766', '2.20', '9.31', '5.76'],
	                   'G': ['0.48', '0', '2', '0', '9', '0', '0.881', '0.179052', '57.0519', '2.35', '9.78', ' 6.06'],
	                   'H': ['-0.4', '-0.5', '4', '79', '10.4', '0.23', '2.025', '-0.01069', '137.1411', '1.80', '9.33', '7.64'],
	                   'I': ['1.38', '-1.8', '2', '93.5', '5.2', '0.186', '1.81', '0.021631', '113.1594', '2.32', '9.76', '6.04'],
	                   'K': ['-1.5', '3', '2', '100', '11.3', '0.219', '2.258', '0.017708', '128.1741', '2.16', '9.06', '9.47'],
	                   'L': ['1.06', '-1.8', '2', '93.5', '4.9', '0.186', '1.931', '0.051672', '113.1594', '2.33', '9.74', '6.04'],
	                   'M': ['0.64', '-1.3', '2', '94.1', '5.7', '0.221', '2.034', '0.002683', '131.1986', '2.13', '9.28', '5.71'],
	                   'N': ['-0.78', '2', '4', '58.7', '11.6', '0.134', '1.655', '0.005392', '114.1039', '2.14', '8.72', '5.43'],
	                   'P': ['0.12', '0', '2', '41.9', '8', '0.131', '1.468', '0.239531', '97.1167', '1.95', '10.64', '6.30'],
	                   'Q': ['-0.85', '0.2', '4', '80.7', '10.5', '0.18', '1.932', '0.049211', '128.1307', '2.17', '9.13', '5.65'],
	                   'R': ['-2.53', '3', '4', '105', '10.5', '0.18', '1.932', '0.049211', '156.1875', '1.82', '8.99', '10.76'],
	                   'S': ['-0.18', '0.3', '4', '29.3', '9.2', '0.062', '1.298', '0.004627', '87.0782', '2.19', '9.21', '5.70'],
	                   'T': ['-0.05', '-0.4', '4', '51.3', '8.6', '0.108', '1.525', '0.003352', '101.1051', '2.09', '9.10', '5.60'],
	                   'V': ['1.08', '-1.5', '2', '71.5', '5.9', '0.14', '1.645', '0.057004', '99.1326', '2.29', '9.74', '6.02'],
	                   'W': ['0.81', '-3.4', '3', '145.5', '5.4', '0.409', '2.663', '0.037977', '186.2132', '2.46', '9.41', '5.88'],
	                   'Y': ['0.26', '-2.3', '3', '117.3', '6.2', '0.298', '2.368', '0.023599', '163.1760', '2.20', '9.21', '5.63']}
chem_df = pd.DataFrame.from_dict(chemical_properties, columns=chemicals_features, orient='index')


# names of species for feature comparison:

species_names = [ \
    'actinomyces_israelii',
    'bacillus_clausii',
    'bacillus_subtilis',
    'escherichia_coli',
    'staph_aureus',
    'strep_mitis',
    'strep_mutans',
    'strep_salivarius',
    'helicobacter_pylori',
    'strep_sanguinis']