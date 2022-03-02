from ComputationalBiology.biore import biore
import re
import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('max_colwidth', 8000)
from sklearn.metrics import confusion_matrix
from ComputationalBiology.biore.biore_macros import PTS, pts_to_run

if __name__ == '__main__':
    if True:
        #compute matches for all regexes
        df = pd.read_csv('./data/human1.tab', sep='\t', skiprows=0)
        print('Total proteins', len(df))
        df['Sequence'] = df['Sequence'].apply(lambda x: x.upper())
        for reg in pts_to_run:
            signal_name = reg.name
            signal_regex = reg.value
            df[signal_name] = df['Sequence'].apply(lambda x: biore.search(signal_regex, x) is not None)
            df[signal_name + '_MATCH'] = df['Sequence'].apply(lambda x: str(biore.search(signal_regex, x)))

        # removing all rows where no label exists (empty string)
        df = df.fillna('')
        print('total rows: ', len(df))
        df['Subcellular location [CC]'].replace('', np.nan, inplace=True)
        df.dropna(subset=['Subcellular location [CC]'], inplace=True)
        print('after removal of rows with no label:', len(df))

        # FILL GROUND TRUTH
        df['isNucleus'] = df['Subcellular location [CC]'].apply(
            # lambda x: re.search(r'\bNucleus\b', x) is not None)
            lambda x: re.search(r'^SUBCELLULAR LOCATION: Nucleus', x) is not None)
        df['isMitochondria'] = df['Subcellular location [CC]'].apply(
            lambda x: re.search(r'^SUBCELLULAR LOCATION: Mitochondrion', x) is not None)
        df['isER'] = df['Subcellular location [CC]'].apply(
            lambda x: re.search(r'^SUBCELLULAR LOCATION: Endoplasmic reticulum', x) is not None)
        df['isHelix'] = df['Helix'].apply(lambda x: x is not '')  # check which helix! alpha/pi

        print('Total nucleus found: ', sum(df['isNucleus'] == True))
        print('Total mitochondria found: ', sum(df['isMitochondria'] == True))
        print('Total ER found: ', sum(df['isER'] == True))

        for reg in pts_to_run:
            conf = confusion_matrix(df['isNucleus'], df[reg.name])
            tn, fp, fn, tp = conf.ravel()
            sensitivity = tp / (tp + fn)
            specificity = tn / (tn + fp)
            print(reg.name, ': ', conf, 'sensitivity: ', sensitivity, 'specificity: ', specificity)
            print('--------------------')

        df.to_pickle('./data/human1.pickle')

    # PRINT LIST OF GENES
    if True:
        df = pd.read_pickle('./data/human1.pickle')
        df.to_csv('./data/test.csv')
        # for reg in pts_to_run:
        #    found = df_exact[df_exact[reg.name] == True]
        #    print(found)

    print('done')
