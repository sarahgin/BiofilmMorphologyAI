from ComputationalBiology.biore import biore
import re
import pandas as pd
pd.set_option('display.max_rows', None)

from ComputationalBiology.biore.biore_macros import PTS, pts_to_run

if __name__ == '__main__':
    # reg = '(NP)*?\K'
    # seq = 'YDREKDREKKDREKDREKK'

    # reg = '\A[PO]{2,5}\Q*?\K\K\K'
    # tr_reg = biore.translate(reg)
    # print(tr_reg)
    # seq = 'ARRKVKKQQQQQQKKKQQQQQKKK'

    # ans = biore.finditer(reg, seq)
    # lst = [x for x in ans]
    # print(lst)

    if True:
        df = pd.read_csv('./data/human1.tab', sep='\t', skiprows=0)
        print('Total proteins', len(df))
        df['Sequence'] = df['Sequence'].apply(lambda x: x.upper())
        for reg in pts_to_run:
            signal_name = reg.name
            signal_regex = reg.value
            df[signal_name] = df['Sequence'].apply(lambda x: biore.search(signal_regex, x) is not None)
        df = df.fillna('')
        df['isNucleus'] = df['Subcellular location [CC]'].apply(lambda x: re.search(r'\bNucleus\b', x) is not None)
        df.to_pickle('./data/human1.pickle')

    #PRINT LIST OF GENES
    if True:
        df_exact = pd.read_pickle('./data/human1.pickle')
        for reg in pts_to_run:
            num_found = len(df_exact[df_exact[reg.name] == True])
            print(reg.name, num_found)

    print('done')
