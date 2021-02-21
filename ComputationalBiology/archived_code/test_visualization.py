print('hello')
import matplotlib.pyplot as plt
import pickle

# from ComputationalBiology.data_analysis.main_parser_features_calc import PREFIX_SUFFIX_DF_FILE

PREFIX_SUFFIX_DF_FILE = '../../data/data_outputs/prefix_suffix_dict_' + 'BS3610' + '.pickle'
print('Loading pickle file: {}...'.format(PREFIX_SUFFIX_DF_FILE))
with open(PREFIX_SUFFIX_DF_FILE, 'rb') as pickle_file:
    prefix_suffix_dict, prefix_count = pickle.load(file=pickle_file)

    # normalize
    for k in prefix_suffix_dict:
        for j in prefix_suffix_dict[k]:
            prefix_suffix_dict[k][j] /= prefix_count[k]
            prefix_suffix_dict[k][j] *= 100

    max_dict = {}
    for k in prefix_suffix_dict:
        d = prefix_suffix_dict[k]
        d_max_val = max(d.values())

        # taking list of car values in v
        v = list(d.values())
        # taking list of car keys in v
        vk = list(d.keys())
        d_max_str = (vk[v.index(max(v))])

        max_dict[k + '-' + d_max_str] = d_max_val

    #sorted_max_dict = {k: v for k1, v1 in sorted(max_dict.items(), key=lambda item: item[1])}

    z = list(range(len(max_dict.values())))
    y = list(max_dict.values())

    fig, ax = plt.subplots()
    ax.scatter(z, y)
    for i, txt in enumerate(max_dict.keys()):
        print(z[i])
        print(y[i])
        print(txt)
        ax.annotate(txt, (z[i], y[i]))

    plt.show()
exit(0)
