#!/usr/bin/env python

"""Tools to convert raw DOLPHOT output to other formats.

May be executed from the command line.
EXTREMELY WORK IN PROGRESS.

"""

import pandas as pd
import re
import vaex

from astropy.io import fits
from astropy.wcs import WCS

__all__ = ['read_colfile', 'ascii_to_vaex', 'select_cols']

# Mappings of DOLPHOT column descriptions to column names
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


def read_colfile(colfile, replace_single=replace_single,
                 replace_global=replace_global):
    """Construct a table of column names for DOLPHOT output.

    Inputs
    ------
    colfile : str, path object or file-like object
        Path to DOLPHOT '.column' file; see pandas.read_csv

    Returns
    -------
    df_col : pandas.DataFrame
        A table of column descriptions and their corresponding names.
    """
    # read column file, splitting into index (column number) and description
    df_col = pd.read_csv(colfile, sep=re.escape('. '), names=['desc'],
                         index_col=0, engine='python')
    # split into "single" and "global" based on having a ', ' before ' ('
    is_single = df_col['desc'].str.split(re.escape(' (')).str[0].\
        str.contains(re.escape(', '))
    df_single = df_col.loc[is_single]
    df_global = df_col.drop(df_single.index)

    # get rid of 'Object' and select first word
    names_global = df_global.desc.str.replace('Object', '').\
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
                                  ).str.replace('^_', '')
    df_col.loc[:, 'names'] = names
    return df_col


def ascii_to_vaex(asciifile, names, usecols=None):
    """Read raw DOLPHOT output photometry into Vaex.

    Inputs
    ------
    asciifile : str, path object or file-like object
        Path to DOLPHOT photometry file; see pandas.read_csv
    names : list-like
        Sequence of column names; see pandas.read_csv
    usecols : list-like or callable, optional
        Subset of columns to be read in; see pandas.read_csv
        If None, assumes 'names' sequence corresponds to first N columns

    Returns
    -------
    ds : vaex.dataframe.DataFrame
        Photometry table.
    """
    usecols = list(range(len(names))) if usecols is None else usecols
    compression = 'gzip' if asciifile.endswith('gz') else 'infer'
    ds = vaex.from_csv(asciifile, copy_index=False, delim_whitespace=True,
                       names=names, usecols=usecols, header=None,
                       na_values=['99.999'], compression=compression)
    return ds


# HERE BE DRAGONS aka non-generalized functions


def select_cols(df_col, regex='(^[X,Y]$)|(^[F,G]Q?[0-9]{2,5}[W,M,N,X,L]P?_)'):
    """Select columns to be used in ascii_to_vaex.

    Only for M33 use right now!!!
    Regex: X + Y + all columns beginning with any ACS/WFC3 filter name
    """
    cols = df_col.loc[df_col.names.str.match(regex), 'names']
    names = cols.tolist()
    usecols = (cols.index - 1).tolist()
    return names, usecols


def conv_pix_wcs(x, y, wcsfile, origin=1, name='all_pix2world'):
    w = WCS(fits.Header.fromtextfile(wcsfile))
    return getattr(w, name)(x, y, origin)


def add_wcs(ds, base, firstcols=['RA', 'DEC', 'X', 'Y']):
    if 'RA' in ds.get_column_names():
        ds.drop('RA', inplace=True, check=False)
    if 'DEC' in ds.get_column_names():
        ds.drop('DEC', inplace=True, check=False)
    wcsfile = f'{base}/{base}_F475W_drc_wcs.txt'
    ra, dec = conv_pix_wcs(*ds.evaluate(['X', 'Y']), wcsfile)
    ds.add_column('RA', ra)
    ds.add_column('DEC', dec)
    names = firstcols + [n for n in ds.get_column_names()
                         if n not in firstcols]
    return ds[names]


def calc_xmed(ds_left, ds_right):
    # get center pixel for overlapping photometry
    return int(round((ds_right.X.min() + ds_left.X.max()) / 2))


def merge_parallel_phot(name, export=True):
    # merge 3 photsec catalogs
    ds1 = vaex.open(f'{name}/{name}_1.phot.hdf5')
    ds2 = vaex.open(f'{name}/{name}_2.phot.hdf5')
    ds3 = vaex.open(f'{name}/{name}_3.phot.hdf5')
    x0, x1 = calc_xmed(ds1, ds2), calc_xmed(ds2, ds3)
    ds1.select(f'X <= {x0}', name='__filter__')
    ds2.select_box(['X'], [[x0, x1]], name='__filter__')
    ds3.select(f'X >= {x1}', name='__filter__')
    ds = vaex.dataframe.DataFrameConcatenated([d.extract() for d in
                                               [ds1, ds2, ds3]])
    ds = add_wcs(ds, name)
    if export:
        ds.export_hdf5(f'{name}/{name}.phot.hdf5')
        print(f'Merged catalogs written to {name}/{name}.phot.hdf5')
    else:
        return ds


if __name__ == '__main__':
    import sys
    import os
    base = sys.argv[1]
    os.system(f'gzip {base}/*.phot')
    for i in range(1, 4):
        photfile = f'{base}/{base}_{i}.phot'
        colfile = photfile + '.columns'
        hdffile = photfile + '.hdf5'
        df_col = read_colfile(colfile)
        names, usecols = select_cols(df_col)
        ascii_to_vaex(photfile + '.gz', names=names,
                      usecols=usecols).export_hdf5(hdffile)
        print(f'Written {hdffile}')
    merge_parallel_phot(base)
