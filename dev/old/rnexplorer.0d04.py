#!/usr/bin/env python3

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.core.cli as corecli
import rotifer.neighborhood.neighborhood as nh
import argparse
import warnings
import numpy as np
warnings.filterwarnings('ignore')
from collections import defaultdict

#TODO: Build class for input, output, filter

__version__ = 0.04
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Convert rneighbors output (table<->gi2operons)',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('file',
                        help = 'Input neighborhood file',
                        action = corecli.action.autoload,
                        nargs = '*',
                        duplicates = False)

    parser.add_argument('-y', '--header',
                        action='append',
                        dest = 'addinfo',
                              help = 'Add more information to gi2operon header. The additional information must be a valid collumn in the table format\n\
(example: -y assembly, -y classification -y seq_type)',
                              default = [])

    parser.add_argument('-of', '--outformat',
                        dest = 'outformat',
                        help = 'Output format (table, gi2ioperon, not implemented)',
                        default = '')

    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rnexplorer',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Convert neighborhood output (table<->gi2operons)'
                                                  ))
    parser.add_argument('--exclude_by_type', '--xt',
                        '-exclude_by_type', '-xt',
                        dest = 'excludebytype',
                        default = [],
                        help = '''Exclude individual features type (usage: KEY:VALUE).
The value could be an unique value (e.g -xt type:CDS) or a file containing one column with values)
example: -xt type:CDS
         -xt seq_type:[file]
         -xt seq_type:plasmid'''
                                 ,
                        action = 'append',
                        type = str)

    parser.add_argument('--include_by_type', '--it',
                        '-include_by_type', '-it',
                        dest = 'includebytype',
                        default = [],
                        help = '''Include features type (usage: KEY:VALUE).
The value could be an unique value (e.g -it type:CDS) or a file containing one column with values).
The key is a valid column in the table.
example: -it type:CDS
         -it seq_type:[file]
         -it seq_type:plasmid''',
                        action = 'append',
                        type = str)
    parser.add_argument('-a', '--annotation',
                        '--a', '-annotation',
                        dest = 'annotation',
                        default = [],
                        help = '''Add annotation column (usage: FILE, COL_KEY:FILE, or COL_KEY:FILE:COL_NAME).
The annotation file must be a two columns file. The first column is the key the second column is the annotation.
There are 3 ways to use this:
Input just the annotation file, It will map the first column against the pid column in the neighborhood table.
example: -a [FILE]
You can pass the column to map and the file containg the annotation.
example: -a pid:[FILE]
         -a nucleotide:[FILE]
You can pass 3 parameters, the first is the column to map, the second the annotation file, and the last parameter is the new column name.
example: -a pid:[FILE]:My_new_col''',
                        action = 'append',
                        type = str)

    args = parser.parse_args()

    return args


def exclude_type(df, exclude_dict):
    for k,v in exclude_dict.items():
        df = df[~(df[k].isin(v))]
    return df

def include_type(df, include_dict):
    for k,v in include_dict.items():
        df = df[df[k].isin(v)]
    return df

def arg2dict(ls_of_str):
    dc = defaultdict(list)
    for select in ls_of_str:
        k,v = select.split(':')
        try:
            a = open(v).read().splitlines()
            dc[k].extend(a)
        except:
            dc[k].append(v)
    return dc

def ann2dict(ls_of_str):
    dc = defaultdict(list)
    for select in ls_of_str:
        splitted = select.split(':')
        if len(splitted) == 1: # map to pid and its a file
            try:
                a = open(splitted[0]).read().splitlines()
                source_ls = []
                target_ls = []
                for line in a:
                    s = line.split('\t')
                    source_ls.append(s[0])
                    target_ls.append(s[1])
                dc['pid'].append(pd.DataFrame({'source':source_ls, 'target':target_ls}))
            except:
                sys.stderr.write('Unable to read {0}\n'.format(splitted[0]))
                sys.exit()

        if len(splitted) == 2: # could be pid:file  or column:file
            try:
                a = open(splitted[1]).read().splitlines()
                source_ls = []
                target_ls = []
                for line in a:
                    s = line.split('\t')
                    source_ls.append(s[0])
                    target_ls.append(s[1])
                dc[splitted[0]].append(pd.DataFrame({'source':source_ls, 'target':target_ls}))
            except:
                sys.stderr.write('Unable to read {0}\n'.format(splitted[1]))
                sys.exit()
        if len(splitted) == 3: # collumn:file:column_name
            # try:
                a = open(splitted[1]).read().splitlines()
                source_ls = []
                target_ls = []
                for line in a:
                    s = line.split('\t')
                    source_ls.append(s[0])
                    target_ls.append(s[1])
                dc[splitted[0]].append((pd.DataFrame({'source':source_ls, 'target':target_ls}), splitted[2])  )
            # except:
            #     sys.stderr.write('Unable to read {0}\n'.format(splitted[1]))
            #     sys.exit()

    return dc


def add_annotation(df, dc):
    new_info_ls = []
    for k in dc.keys(): # k is the column
        for ann_df in dc[k]: # each df
            if len(ann_df) ==1:
                try: # try to merge
                    df = df.merge(ann_df, how = 'left', left_on = k, right_on = 'source')
                    df.drop(columns= ['source'], inplace = True)
                    df.columns = (list(df.columns[:-1])+[str(len(df.columns))])
                    df = df.fillna('.')
                    new_info_ls.append(str(len(df.columns)))
                except: pass
            if len(ann_df) == 2:
                try: # try to merge
                    ann, col_ann = ann_df

                    df = df.merge(ann, how = 'left', left_on = k, right_on = 'source')
                    df.drop(columns= ['source'], inplace = True)
                    df.columns = (list(df.columns[:-1])+[col_ann])
                    df = df.fillna('.')
                    new_info_ls.append(col_ann)
                except: pass

    return (new_info_ls,df)
#     df[str(len(df.collumns))] #left join


class filter:
    def __init__(self):
        pass

if __name__ == '__main__':
    args = parse_cli()
    outformat = args.outformat
    if outformat:
        ofl = args.outformat.lower()
        if ofl in ['gi2operon','gi2operons']:
            outformat = 'gi2operon'
        elif ofl in ['table', 'tables']:
            outformat = 'table'
        else:
            sys.stderr.write('Unknown option {0}\n'.format(args.outformat))
            sys.stderr.write('Select a valid output format')
            sys.exit()
    ipt = nh.ipt2table(args.file)
    df = ipt.df
    new_info_ls = ''

    # Additional information

    other_info = args.addinfo

    # Include/Exclude Block

    exclude_dict = defaultdict(list)
    include_dict = defaultdict(list)
    excludebytype = args.excludebytype
    includebytype = args.includebytype
    if excludebytype:
        exclude_dict = arg2dict(excludebytype)
    if includebytype:
        include_dict = arg2dict(includebytype)

    # Annotation Block

    annotation_dict = defaultdict(list)
    annotation = args.annotation

    if annotation:
        annotation_dict = ann2dict(annotation)
    # Select output format
    if outformat: pass
    else:
        outformat = 'gi2operon' if ipt.input_format_type == 'gi2operon' else 'table'

    # Apply exclude filter
    if exclude_dict and not include_dict:
        df = exclude_type(df, exclude_dict)

    # Apply include filter
    if include_dict:
         df = include_type(df, include_dict)

    # Apply annotation dict
    if annotation_dict:
        new_info_ls, df = add_annotation(df, annotation_dict)
        # Merge
    if not df.empty:
        nh.printer(df, outformat, other_info, new_info_ls)
