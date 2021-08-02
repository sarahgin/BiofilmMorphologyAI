import logomaker
import numpy as np
import pandas as pd

from ComputationalBiology.biore.biore_macros import AA


def create_logo(items):

    df = pd.DataFrame(items, columns=['Seq'])

    for p in range(21):
        df['position_' + str(p)] = df['Seq'].apply(lambda x: x[p])

    df_logo = pd.DataFrame(np.zeros((21, 20)), columns=[e.value for e in AA])

    for p in range(21):
        for g in [e.value for e in AA]:
            df_logo.loc[p, g] = (sum(df['position_' + str(p)] == g))

    seq_num = df_logo.sum(axis=1).iloc[0]
    df_logo = df_logo.applymap(lambda i: i / seq_num)

    # create Logo object
    ww_logo = logomaker.Logo(df_logo,
                             font_name='Stencil Std',
                             color_scheme='NajafabadiEtAl2017',
                             vpad=.1,
                             width=.8)

    # style using Logo methods
    ww_logo.style_xticks(anchor=0, spacing=5, rotation=45)
    ww_logo.highlight_position(p=4, color='gold', alpha=.5)
    ww_logo.highlight_position(p=26, color='gold', alpha=.5)

    # style using Axes methods
    #ww_logo.ax.set_ylabel('information (bits)')
    #ww_logo.ax.set_xlim([-1, len(ww_df)])
    return ww_logo
