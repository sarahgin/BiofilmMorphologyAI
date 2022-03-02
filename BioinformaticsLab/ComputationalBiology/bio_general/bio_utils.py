from itertools import product

from BioinformaticsLab.ComputationalBiology.bio_general.bio_macros import ValidAlphabet, alphabets_dict


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


def merge_add_dicts(dict1, dict2):
    merged_dict = {key: dict1.get(key, 0) + dict2.get(key, 0)
                   for key in set(dict1) | set(dict2)}
    return merged_dict
