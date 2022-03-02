import pandas as pd
import numpy as np

# load annotations:
df_func = pd.read_csv('C:/Users/user/Downloads/geneName.txt')  # attached - this is the column locus tag from geneCatrgories.csv
print(len(df_func))
print(df_func.head())

unique_annotation = df_func['locus tag'].unique()
print('num unique annotations:', len(unique_annotation))

# load our pickle:
species_name = 'BS3610'
FEATURES_DF_FILE = '../../data/data_outputs/features_' + species_name + '.pickle'  # attached
species_df = pd.read_pickle(FEATURES_DF_FILE)
print(species_df.columns)
print(len(species_df))

unique_gb = species_df['GENE_NAME'].unique()
print('num unique genes in genBank:', len(unique_gb))

# find intersection:
intersect = np.intersect1d(unique_annotation, unique_gb)
print(len(intersect))


from collections import Counter
dict_counts = Counter(species_df['GENE_NAME'])
print(dict_counts)

# find repeating genes:
repeating_genes = {}
num_repeats = 0
for k in dict_counts:
    if dict_counts[k] > 1 and k != '':
        repeating_genes[k] = dict_counts[k]
        num_repeats += dict_counts[k]

print('repeating_genes:', repeating_genes)
print('num_repeats:', num_repeats)


