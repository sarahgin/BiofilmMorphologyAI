import os

from ComputationalBiology.file_utils import genbank_parser

if __name__ == '__main__':
    genbank_file = '../../data/GeneticData/BS3610.gb'
    assert(os.path.exists(genbank_file))  # making sure that the path is valid
    record_gb = genbank_parser.read_genbank_file(genbank_file)
    genbank_parser.init_species(record_gb)

    print(genbank_parser.get_stats(record_gb))
    print(record_gb)

# TODO: operons?
# TODO: assert that locus_tag is unique
#
# # for cds: use codon_start, translation_table

