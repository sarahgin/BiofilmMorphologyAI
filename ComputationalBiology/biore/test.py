from ComputationalBiology.biore import biore

if __name__ == '__main__':
    my_regex = '(NP)*?\K'
    protein1 = 'YDREKDREKKDREKDREKK'

    #ans = biore.search(my_regex, protein1)
    #print(ans)

    ans = biore.finditer(my_regex, protein1)
    lst = [x for x in ans]
    print(lst)