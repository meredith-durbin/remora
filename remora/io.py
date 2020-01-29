#!/usr/bin/env python

import pandas as pd
import re

replace_single = {
    '((Total)|(Measured)) counts'       : 'COUNT',
    '((Total)|(Measured)) sky level'    : 'SKY',
    'Normalized count rate uncertainty' : 'RATERR',
    'Normalized count rate'             : 'RATE',
    'Instrumental VEGAMAG magnitude'    : 'VEGA',
    'Transformed UBVRI magnitude'       : 'TRANS',
    'Magnitude uncertainty'             : 'ERR',
    'Photometry quality flag'           : 'FLAG',
}

replace_global = {
    'Signal-to-noise' : 'SNR',
    '(ness)|(ing)'    : '',
}

def get_colnames(colfile, replace_single=replace_single,
                 replace_global=replace_global):
    """Construct a table of column names for DOLPHOT output, with indices
    corresponding to the column number in the DOLPHOT output file.

    Inputs
    ------
    colfile : str, path object or file-like object
        Path to DOLPHOT '.column' file; see pandas.read_csv

    Returns
    -------
    df_col : DataFrame
        A table of column descriptions and their corresponding names.
    """
    # read column file, splitting into index (column number) and description
    df_col = pd.read_csv(colfile, sep='\.\ ', names=['desc'],
                         index_col=0, engine='python')
    # split into "single" and "global" based on having a ', ' before ' ('
    is_single = df_col['desc'].str.split(re.escape(' (')).str[0].\
                               str.contains(re.escape(', '))
    df_single = df_col.loc[is_single]
    df_global = df_col.drop(df_single.index)
    
    names_global = df_global.desc.str.replace('Object','').\
                                  str.strip().str.split().str[0]
    
    single_split = df_single['desc'].str.split(re.escape(', '))
    
    names_single = single_split.str[0]
    for k, v in replace_single.items():
        names_single = names_single.replace(re.compile(k), v)
    
    names_concat = pd.concat([names_global, names_single]
                            ).reindex_like(df_col)
    for k, v in replace_global.items():
        names_concat = names_concat.replace(re.compile(k), v)
        
    prefix_global = pd.Series('', index=names_global.index)
    prefix_single = single_split.str[1].str.split().str[0].\
                                 str.replace(re.escape('.chip'), '_chip')
    prefix_concat = pd.concat([prefix_global, prefix_single]
                             ).reindex_like(df_col)
    names = prefix_concat.str.cat(names_concat.str.upper(), sep='_'
                                 ).str.replace('^_','')
    df_col.loc[:, 'names'] = names
    return df_col
