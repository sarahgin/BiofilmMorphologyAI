from ComputationalBiology.biore.biore_macros import AA_GROUP, group_to_aa_dict, aa_group_list, aa_to_group_dict
import re


def translate(original_str: str):
    res = ''
    i = 0
    while i < len(original_str):
        c = original_str[i]
        if c == "\\":
            res += original_str[i + 1]
            i += 1
        elif c in aa_group_list:
            res += group_into_aa_or(c)
        else:
            res += c
        i += 1
    return res


def group_into_aa_or(group_letter: str):
    # TODO: handle case of invalid group letter
    relevant_aa = group_to_aa_dict[group_letter]
    res = '('
    for i in range(len(relevant_aa) - 1):
        res += relevant_aa[i] + '|'
    res += relevant_aa[i + 1] + ')'
    return res


# simple translations single-letter aa into amino acid group letter
def aa_into_group(seq):
    res = ''
    for c in seq:
        res += aa_to_group_dict[c]
    return res


def group_into_aa(seq):
    res = ''
    for c in seq:
        res += group_to_aa_dict[c]
    return res
