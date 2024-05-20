from BioinformaticsLab.ComputationalBiology.biore.biore_macros import lehninger_group_to_aa_dict, \
    lehninger_aa_group_list, LEHNINGER_AA_GROUP


def translate_biore_into_re(original_str: str):
    '''
    Given string composed of either a lehninger group or an amino acid,
    the function returns a translation.
    e.g., 'P' will be translated into (K|R|H), '\P' will be translated into P (Proline)
    '''
    res = ''
    i = 0
    while i < len(original_str):
        c = original_str[i]
        if c == "\\":
            res += original_str[i + 1]
            i += 1
        elif c in lehninger_aa_group_list:
            res += lehninger_group_into_aa_or(c)
        else:
            res += c
        i += 1
    return res


def lehninger_group_into_aa_or(group_letter: str):
    '''
    Translates a group symbol to a regex of actual amino acids separated by or.
    e.g., 'P' (positive group) will be translated into (K|R|H)
    '''
    # TODO: handle case of invalid group letter
    relevant_aa = lehninger_group_to_aa_dict[group_letter]
    res = '('
    for i in range(len(relevant_aa) - 1):
        res += relevant_aa[i] + '|'
    res += relevant_aa[i + 1] + ')'
    return res


if __name__ == '__main__':
    res = translate_biore_into_re('\P\P\P.{1,3}\L\L\L')
    print(res)

