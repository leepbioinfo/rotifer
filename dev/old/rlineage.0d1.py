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

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Return lineage from a NCBI lineage column')
    parser.add_argument('dataframe',

                        help = 'Input a table',
                        nargs = '*',
                        action = corecli.action.autoload, duplicates = False)

    parser.add_argument('--outformat', '-o', '--of', '-of', dest = 'outformat',
                        help = 'Output format (full, partial)',
                        default = 'full')

    parser.add_argument('--column', '-c', '--c',
                        dest = 'column',
                        default = 'classification',
                        help = 'Input column name containing NCBI lineage information (Default = classification)')

    parser.add_argument('-y', '--y',
                        '-header', '--header',
                        dest = 'header',
                        action = 'append',
                        help = 'Output columns header (default: Classification)',
                        default = ['all'],
                        type = str)

    parser.add_argument('-t', '--taxonomy',
                        '--t', '-taxonomy',
                        dest = 'taxonomy',
                        action = 'append',
                        help = 'Add taxonomic group in the output (Default: superkingdom, phylum, class) [Options: superkingdom, phylum, class, order, family, genus, species]',
                        default = []
                        )

    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rlineage',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Return lineage'))

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



