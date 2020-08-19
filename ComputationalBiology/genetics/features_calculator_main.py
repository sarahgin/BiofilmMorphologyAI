import os

from BiofilmMorphologyAI.ComputationalBiology.biologyutils.WholeGenomeSequence import SpeciesGenomeSequence
from BiofilmMorphologyAI.ComputationalBiology.fileutils.gene_bank_parser import \
    filter_tuples, genbank_to_tuples, read_genbank_file
from BiofilmMorphologyAI.ComputationalBiology.genetics.features_calculator import get_tuple_features

if __name__ == '__main__':
    # open genes file
    gene_bank_file = '../../data/GeneticData/BS3610.gb'
    assert(os.path.exists(gene_bank_file))  # making sure that the path is valid
    record_gb = read_genbank_file(gene_bank_file)
    lst, all_types = genbank_to_tuples(record_gb)

    # list of genes that code to proteins, and their matching protein
    lst_CDS = filter_tuples(lst, product_type='CDS')

    # TODO: handle minus strand
    genome = SpeciesGenomeSequence(record_gb.seq)
    # # 1st gene
    # start1 = lst_CDS[0].gene.location.start.position
    # end1 = lst_CDS[0].gene.location.end.position
    # strand1 = lst_CDS[0].gene.location.strand
    # gene_sequence1 = genome.get_reference_seq(start1, end1)
    # aa_num = len(lst_CDS[0].gene_product.qualifiers['translation'][0])

    # gene_seq from the genome
    # if strand == -1:
    # gene_seq = rc(gene_seq)
    # gene_seq[0:3] == 'ATG' and gene_seq[-3:] in ['TAG', 'TGA', 'TAA]
    # and (len(gene_seq) - 3) == num_aa * 3

    for i in range(len(lst_CDS)):
        start1 = lst_CDS[i].gene.location.start.position
        end1 = lst_CDS[i].gene.location.end.position
        strand1 = lst_CDS[i].gene.location.strand
        gene_seq = genome.get_coding_sequence(start1, end1, strand1)
        print(i)
        protein_seq = lst_CDS[i].gene_product.qualifiers['translation'][0]
        # maybe operon?

        if not(gene_seq[0:3] == 'ATG' and gene_seq[-3:] in ['TAG', 'TGA', 'TAA']):
            print((len(gene_seq) - 3), len(protein_seq) * 3)
            print(lst_CDS[i])


        # if lst_CDS[i].gene.location.strand == -1:
        #     print(i)
        #     start2 = lst_CDS[i].gene.location.start.position
        #     end2 = lst_CDS[i].gene.location.end.position
        #     gene_sequence2 = genome.get_reference_seq(start2, end2)
        #     print(gene_sequence2[-3:])
        #     lst_CDS[i].qualifiers['translation']
        #     print()

    # for a single gene+protein:
    dict_features = get_tuple_features(lst_CDS[0], SpeciesGenomeSequence(record_gb.seq))
    # compute GC content for 1st gene
    print(dict_features)



    # TODO: compute to all, save in dataFrame of: rows = genes, cols = features

    # # compute the features per gene
    # lst_features = [] # TODO: convert to dataFrame
    #
    # # compute avg. over the genes



# if __name__ == '__main__':
#     gene_bank_file = '../../data/GeneticData/BS3610.gb'
#     assert(os.path.exists(gene_bank_file))  # making sure that the path is valid
#     record_gb = read_gene_bank_file(gene_bank_file)
#     # print(record_gb.annotations)
#     # print(type(record_gb.annotations))
#
#     features = record_gb.features
#     print(record_gb.features)
#     print(type(record_gb.features))
