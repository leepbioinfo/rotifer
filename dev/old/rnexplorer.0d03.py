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

__version__ = 0.03
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


class filter:
    def __init__(self):
        pass

if __name__ == '__main__':
    args = parse_cli()
    # dates_dict = defaultdict(list)
    # for key, date in cur:
    #         dates_dict[key].append(date)
    # to_open = open(args.file)
    # Check option exit if invalid
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

    other_info = args.addinfo

    exclude_dict = defaultdict(list)
    include_dict = defaultdict(list)
    excludebytype = args.excludebytype
    includebytype = args.includebytype
    if excludebytype:
        for select in excludebytype:
            k,v = select.split(':')
            try:
                a = open(v).read().splitlines()
                exclude_dict[k].extend(a)
            except:
                exclude_dict[k].append(v)
    if includebytype:
        for select in includebytype:
            k,v = select.split(':')
            try:
                a = open(v).read().splitlines()
                include_dict[k].extend(a)
            except:
                include_dict[k].append(v)

    if outformat:
        pass
    else:
        outformat = 'gi2operon' if ipt.input_format_type == 'gi2operon' else 'table'

    if exclude_dict and not include_dict:
        df = exclude_type(df, exclude_dict)
    if include_dict:
         df = include_type(df, include_dict)

    nh.printer(df, outformat, other_info)
