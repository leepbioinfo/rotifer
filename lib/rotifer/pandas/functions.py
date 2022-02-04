#!/usr/bin/env python3

import sys
import pandas as pd
from collections import OrderedDict

def optimize_memory_usage(df):
    '''
    Optimize memory usage of a Pandas DataFrame.
    
    Usage:
      from rotifer.pandas.functions import optimize_df
      df = optimize_memory_usage(df)

    Returns:
      A memory optimized version of the input Pandas
      DataFrame with the exact same content as the input.

    Parameters:
      df : a Pandas DataFrame
    '''

    # Make a copy of the input DataFrame
    df = df.copy()

    # int optimization
    converted_int = df.select_dtypes(include=['int'])
    converted_int = converted_int.apply(pd.to_numeric, downcast = "integer")

    # Float
    converted_float = df.select_dtypes(include=['float'])
    converted_float = converted_float.apply(pd.to_numeric,downcast='float')

    # Object
    df_obj = df.select_dtypes(include=['object'])
    converted_obj = pd.DataFrame()

    for col in df_obj.columns:
        num_unique_values = len(df_obj[col].unique())
        num_total_values = len(df_obj[col])
        if num_unique_values / num_total_values <= 0.5:
            converted_obj.loc[:,col] = df_obj[col].astype('category')
        else:
            converted_obj.loc[:,col] = df_obj[col]

    # Adding converted objects to a better df
    df[converted_int.columns] = converted_int
    df[converted_float.columns] = converted_float
    df[converted_obj.columns] = converted_obj
    return df

def print_everything(max_rows=10000000000, max_columns=1000000000, max_colwidth=1000000000, width=100000, verbose=False):
    """
    A function to set pandas display options to print entire
    DataFrames instead of truncating output to a few columns
    and rows.

    Usage:
      import rotifer.pandas.functions as rpf
      rpf.print_everything()

    Parameters:
      max_rows     : maximum number of DataFrame rows to print
      max_columns  : maximum number of DataFrame columns to show
      max_colwidth : maximum length of columns
                     All columns wider than this will be truncated
      width        : maximum length for rows
      verbose      : warn on STDERR when changing options
    
    WARNING:
      Very large Pandas DataFrames could take too long to print!!!
    """

    pd.options.display.max_rows = max_rows
    pd.options.display.max_columns = max_columns
    pd.options.display.max_colwidth = max_colwidth
    pd.options.display.width = width

    if verbose:
        print('Pandas options have been reset!', file=sys.stderr)

def show_display_options():
    '''
    List all Pandas display options that control printed output
    formatting (see pandas.options.display).
    '''
    return pd.DataFrame([ (x,getattr(pd.options.display,x)) for x in dir(pd.options.display) ], columns=('option','value'))

def filter_nonoverlapping_regions(df, reference=['sequence'], start='estart', end='eend', maximum_overlap=0.1, criteria = OrderedDict({'evalue':True, 'region_length':False})):
    '''
    This function selects the best set of non-overlapping intervals.

    The input Pandas DataFrame must conatins columns corresponding
    to unidimensional interval coordinates (start and end columns)
    and at least one reference column.

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    reference : (list of) strings, default ['sequence']
        Name of the column storing sequence identifiers
    maximum_overlap : integer, default 0.1 (10% of the shortest)
        Maximum number of overlapping residues to ignore.
        If maximum_overlap < 1, the parameter is treated as a percentage 
    start : string, default estart
        Name of the column with the first coordinate of each interval
    end : string, default eend
        Name of the column with the last coordinate of each interval
    criteria : ordered dictionary
        Criteria for best interval selection.

        Keys in this dictionary must match column names in the input
        Pandas DataFrame. Values should be booleans (True or False).

        'region_length' is an exception: this column is automatically
        calculated from the start and end coordinates of each region.

        A True value for a given column will give precedence for rows
        with the smallest value in that column. A False value
        will select rows with the largest values.

    Examples
    --------
    HMMER results, parsed by hmmer2table, using the environment limits:

    >>> df = pd.read_csv("results.hmmer2table.tsv", sep="\t")
    >>> dfc = filter_overlapping_regions(df)

    HH-suite results based on the largest probability and/or region length:

    >>> df = pd.read_csv("results.hhsuite2table.tsv", sep="\t")
    >>> c = {'probability':False, 'region_length':False}
    >>> dfc = filter_overlapping_regions(df, start='qstart', end='qend', criteria=c)
    '''
    import numpy as np

    # Process arguments and select columns for analysis
    if not isinstance(criteria,OrderedDict):
        criteria = OrderedDict(criteria)
    if not isinstance(reference,list):
        reference = [ reference ]

    # Copy and sort the dataframe
    # Sort row priority: best rows must appear first
    crit = list(criteria.keys())
    clean = df.copy()
    if 'region_length' in criteria:
        if 'region_length' in clean.columns:
            clean.rename({'region_length':'_original_region_length'}, axis=1, inplace=True)
        clean['region_length'] = clean[end] - clean[start] + 1
    clean.sort_values(reference + crit, ascending=len(reference)*[True] + [ criteria[x] for x in crit ], inplace=True)

    # Remove non-optimal layers
    keep = []
    while len(clean) > len(keep):
        s = clean.loc[clean.index.drop(keep)].drop_duplicates(reference)
        keep.extend(s.index.tolist())
        s = clean.loc[clean.index.drop(keep)].reset_index().rename({'index':'_oid'}, axis=1).merge(s, on=reference, how='inner')
        overlap  = np.where(s[f'{end}_x'] < s[f'{end}_y'], s[f'{end}_x'], s[f'{end}_y'])
        overlap -= np.where(s[f'{start}_x'] > s[f'{start}_y'], s[f'{start}_x'], s[f'{start}_y'])
        overlap += 1
        maxoverlap = maximum_overlap
        if maxoverlap < 1:
            maxoverlap = np.floor(maxoverlap * np.where(s.region_length_x < s.region_length_y, s.region_length_x, s.region_length_y))
        s = s[overlap > maxoverlap]._oid.unique().tolist()
        clean.drop(s, inplace=True)

    # Remove region_length and return
    if 'region_length' in criteria:
        clean.drop(['region_length'], axis=1, inplace=True)
        if '_original_region_length' in clean.columns:
            clean.rename({'_original_region_length':'region_length'}, axis=1, inplace=True)
    clean.sort_values(reference + [start,end], ascending=(len(reference)+2)*[True], inplace=True)
    return clean
