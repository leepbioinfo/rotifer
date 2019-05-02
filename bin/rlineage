#!/usr/bin/env python3

import sys
sys.setrecursionlimit(10000)
sys.path.insert(0,'/home/kaihami/mymodules/')
from rotifer import taxonomy as taxon
import rotifer.neighborhood.neighborhood as nh
import pandas as pd
import numpy as np
import argparse
import rotifer.core.cli as corecli

__version__ = '0.1'
__authors__ = 'Gilberto Kaihami'
__rdoc__='''
DESCRIPTION:
Get lineage
'''

def parse_cli():
    parser = corecli.parser(description = 'get lineage')
    parser.add(dest = 'dataframe',
                        helper = 'input a table',
                        nargs = '*',
                        action = corecli.action.autoload, duplicates = false)

    parser.add('--column', '-c',
               dest = 'column',
               arg_type = str,
                        default = 'classification',
                        helper = 'input column name containing ncbi lineage information (default = classification)')

    parser.add('--header', '-y',
                        dest = 'header',
                        action = 'append',
                        helper = 'output columns header (default: classification)',
                        default = ['all'],
                        arg_type = str)

    parser.add('-t', '--taxonomy',
                        dest = 'taxonomy',
                        action = 'append',
                        helper = 'add taxonomic group in the output (default: superkingdom, phylum, class) [options: superkingdom, phylum, class, order, family, genus, species]',
                        default = []
                        )

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()
    header = args.header

    if args.taxonomy:
        taxonomy = [x.lower() for x in args.taxonomy]
    else:
        taxonomy = ['superkingdom','phylum', 'class']

    df = nh.ipt2table(args.dataframe).df

#    print(header)
    df_header = list(df.columns)
    if len(header) == 1:
        header = df_header + ['lineage']
    else:
        header = [x for x in header if x != 'all']
#        print(header)

    df = df.reset_index()

    tax_df = taxon.taxonomy.df2tax(df).rlineage(taxonomy)
    tax_df = tax_df.reset_index()

    s2 = df.merge(tax_df, how = 'left', left_on = 'index',
                  right_on ='idx')

    s2 = s2[header]
    s2 = s2.fillna('.')
    s2 = s2.applymap(str)
    print('\t'.join(s2.columns))
    for i, row in s2.iterrows():
        print('\t'.join(row.values))



