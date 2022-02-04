#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from collections import OrderedDict

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
