from ComputationalBiology.biore import biore
import re
import pandas as pd
pd.set_option('display.max_rows', None)

from ComputationalBiology.biore.biore_macros import PTS_EXACT, protein_targeting_signals_exact_dict, PTS_GENERAL, \
    protein_targeting_signals_dict

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
        df['Sequence'] = df['Sequence'].apply(lambda x: x.upper())
        for exact_reg in PTS_EXACT:
            signal_name = exact_reg.name
            signal_regex = protein_targeting_signals_exact_dict[exact_reg]
            df[signal_name] = df['Sequence'].apply(lambda x: biore.search(signal_regex, x) is not None)
        pd.to_pickle(df, './data/human1_exact.pickle')

        for general_reg in PTS_GENERAL:
            signal_name = general_reg.name
            signal_regex = protein_targeting_signals_dict[general_reg]
            df[signal_name] = df['Sequence'].apply(lambda x: biore.search(signal_regex, x) is not None)

        df = df.fillna('')
        df['isNucleus'] = df['Subcellular location [CC]'].apply(lambda x: re.search(r'\bNucleus\b', x) is not None)
        pd.to_pickle(df, './data/human1_general.pickle')

    #PRINT LIST OF GENES
    if True:

        df_exact = pd.read_pickle('./data/human1_exact.pickle')
        for exact_reg in PTS_EXACT:
            found = len(df_exact[df_exact[exact_reg.name] == True])
            print(exact_reg.name, found)
        print('------------------------------------------------')
        df_general = pd.read_pickle('./data/human1_general.pickle')
        for general_reg in PTS_GENERAL:
            found = len(df_general[df_general[general_reg.name] == True])
            print(general_reg.name, found)

    print('done')
