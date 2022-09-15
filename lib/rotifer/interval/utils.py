#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from collections import OrderedDict
_criteria = OrderedDict({'sum_probs':False, 'probability':False, 'score':False, 'region_length':False})

def architecture(df, reference=['query'], collapse='hit', start='qstart', end='qend', maximum_overlap=0.1, criteria=_criteria):
    '''
    This function builds a compact representation of all non-overlapping
    regions for a given reference column.

    It is often used to produce a compact protein domain architecture table

    The input Pandas DataFrame must contain columns corresponding
    to unidimensional interval coordinates (start and end columns)
    and at least one reference column (e.g. pid or sequence).

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    reference : (list of) strings, default ['query']
        Name of the column storing sequence identifiers
    collapse : (list of) strings, default 'hits'
        Name of the column storing region identifiers
    start : string, default estart
        Name of the column with the first coordinate of each interval
    end : string, default eend
        Name of the column with the last coordinate of each interval
    maximum_overlap : integer, default 0.1 (10% of the shortest)
        Maximum number of overlapping residues to ignore.
        If maximum_overlap < 1, the parameter is treated as a percentage 
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
    HH-suite results based on the largest probability and/or region length:

    >>> from rotifer.interval.utils import architecture
    >>> df = pd.read_csv("results.hhsuite2table.tsv", sep="\t")
    >>> dfc = architecture(df, start='qstart', end='qend', criteria=c)
    '''
    arch = filter_nonoverlapping_regions(df, reference=reference, start=start, end=end, maximum_overlap=maximum_overlap, criteria=criteria)
    arch.sort_values([*reference, start, end], ascending=True, inplace=True)
    arch = arch.groupby(reference).agg(architecture=(collapse,lambda x: "+".join([ str(y) for y in x ]))).reset_index()
    return arch

def filter_nonoverlapping_regions(df, reference=['sequence'], start='estart', end='eend', maximum_overlap=0.1, criteria = OrderedDict({'evalue':True, 'region_length':False})):
    '''
    This function selects the best set of non-overlapping intervals.

    The input Pandas DataFrame must contain columns corresponding
    to unidimensional interval coordinates (start and end columns)
    and at least one reference column (e.g. pid or sequence).

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    reference : (list of) strings, default ['sequence']
        Name of the column storing sequence identifiers
    start : string, default estart
        Name of the column with the first coordinate of each interval
    end : string, default eend
        Name of the column with the last coordinate of each interval
    maximum_overlap : integer, default 0.1 (10% of the shortest)
        Maximum number of overlapping residues to ignore.
        If maximum_overlap < 1, the parameter is treated as a percentage 
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

    # Process arguments and select columns for analysis
    if not isinstance(criteria,OrderedDict):
        criteria = OrderedDict(criteria)
    if not isinstance(reference,list):
        reference = [ reference ]

    # Copy and sort the dataframe
    # Sort row priority: best rows must appear first
    crit = list(criteria.keys())
    clean = df.copy()
    clean['__rl'] = clean[end] - clean[start] + 1 # Added by Gian to quickly fix the code
    if 'region_length' in clean.columns:
        remove_region_length = False
    else:
        remove_region_length = True
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
            maxoverlap = np.floor(maxoverlap * np.where(s.__rl_x < s.__rl_y, s.__rl_x, s.__rl_y))
        s = s[overlap > maxoverlap]._oid.unique().tolist()
        clean.drop(s, inplace=True)

    # Remove region_length and return
    clean.drop(['__rl'], axis=1, inplace=True) # added by Gian to quickly fix the code
    if remove_region_length:
        clean.drop(['region_length'], axis=1, inplace=True)
    clean.sort_values(reference + [start,end], ascending=(len(reference)+2)*[True], inplace=True)
    return clean
