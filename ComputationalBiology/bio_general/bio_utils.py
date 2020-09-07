from itertools import product

from ComputationalBiology.bio_general.bio_macros import ValidAlphabet, alphabets_dict


def kmers_generator(k: int, alphabet: ValidAlphabet):
    """
    This is a generator, generates all DNA k-mers
    :param alphabet:
    :return: the next k-mer (using yield)
    :param k:
    """
    gen = product(alphabets_dict[alphabet], repeat=k)
    for x in gen:
        yield ''.join(x)


def get_all_kmers(k, alphabet):
    all_kmers = []
    for kmer in kmers_generator(k, alphabet):
        all_kmers.append(kmer)
    return all_kmers