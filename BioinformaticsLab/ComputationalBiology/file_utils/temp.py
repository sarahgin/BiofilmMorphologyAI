from Bio import SeqIO

if __name__ == '__main__':
    gene_bank_file =''
    with open(gene_bank_file, 'r') as input_handle:

        gen = SeqIO.parse(input_handle, 'genbank')
        record_gb = next(gen)  # content of 1st record
        print(record_gb)