from functools import reduce
import pandas as pd
from ch_base import *


def ch_column_value_replace(db, ch_column, ch_replacer):
    db[ch_column] = db[ch_column].str.replace(
        ch_replacer[0], ch_replacer[1])
    return db


def ch_partial_melt(df, melt_column, change_columns, save_columns='all', suffixdict=None, suffix_delimeter='.', method='merge'):
    if type(change_columns) != list:
        raise TypeError('change_columns must be a list')

    if save_columns != 'all' and type(save_columns) != list:
        raise TypeError('save_columns must be a string "all", list or None')
    elif type(save_columns) == list:
        wcolumns = save_columns
    else:
        wcolumns = list(set(df.columns)-set(change_columns)-set([melt_column]))
    
    if suffixdict is None:
        melt_col_values = pd.unique(df[melt_column])
        suffixdict = {key:key for key in melt_col_values}
    
    df_dict = {key:df[df[melt_column]==key][wcolumns+change_columns] for key in suffixdict}
    for key in df_dict:
        df_dict[key].columns = wcolumns + ['{}{}{}'.format(change_column_name, suffix_delimeter, suffixdict[key]) for change_column_name in change_columns]

    if method == 'merge':
        return reduce(lambda left, right: pd.merge(left, right, on=wcolumns, how='outer'),
                      [df_dict[key] for key in df_dict])

    elif method == 'append':
        return reduce(lambda left, right: left.append(right, ignore_index=True),
                      [df_dict[key] for key in df_dict])