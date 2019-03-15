#!/usr/bin/env python3

import os
import subprocess
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.core.cli as corecli
import urllib.request
from Bio import SwissProt
import pandas as pd
import argparse
from termcolor import colored
import argcomplete
from multiprocessing import Pool, Process

__version__ = '0.04'
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser()

    parser.add(rconfig = ':cli.acc')

    parser.add(short_arg = '-o',
                long_arg = '--output',
                helper = 'Output redundant (all) or non-redundant accessions (nr) [Default = all]',
                default = 'all')
    # parser.add(long_arg = '--version',
    #             action = 'version',
    #                     version = corecli.version(program = 'uniprot2gi',
    #                                               version = __version__,
    #                                               authors = __authors__,
    #                                               description = 'Convert UniprotID to NCBI accession')
    #                     )

    parser.add(long_arg = '-nr',
                        helper = 'Return non-redundant accessions (same effect as -o nr)',
                        action = 'store_true')

    parser.add(short_arg = '-t',
               long_arg = '--threads',
                default = 5,
                arg_type = int,
                helper = 'Number of threads (Default: 5)')

    argcomplete.autocomplete(parser)

    args = parser.parse_args()

    return args

def convert(ls, of = 'all'):
    df = pd.DataFrame()
    for ele in ls:
        uniprot_id = ele
        modified = []
        originals = []
        try:
            url = 'http://www.uniprot.org/uniprot/'+uniprot_id+'.txt'
            s = urllib.request.urlopen(url)
            contents = s.read().decode('utf-8')
            contents_format = contents.split('\n')
            for ele in contents_format:
                if 'RefSeq' in ele:
                    found = ele.split(';')[1:]
                    found = [x.strip() for x in found]
                    found2 = []
                    for ele in found:
                        if ele.endswith('.'):
                            found2.append(ele[:-1])
                        else:
                            found2.append(ele)
            modified.extend(found2)
            originals.extend([uniprot_id]*len(found2))
            df = pd.DataFrame({'original': originals,
                                'modified': modified})
            if of == 'all':
                for e in df[['original', 'modified']].values:
                    print('\t'.join(e))
            if of == 'nr':
                df.drop_duplicates('original', inplace = True)
                for e in df[['original', 'modified']].values:
                    print('\t'.join(e))
                    sys.stdout.flush()

        except:
            print('\t'.join([ele,ele]))
            sys.stdout.flush()
    return df


def version():
    s = '''{a}
{program}
{description}
{version}
{authors}
{a}
'''.format(a = '#'*20,
            program = colored('uniprot2gi:','red',attrs=['bold']),
            description = colored('Convert UniprotID to NCBI accession', 'green'),
           version = 'Current version: {0}'.format(colored(__version__, 'cyan')),
           authors = 'Authors: {}'.format(colored(__authors__, 'cyan'))
            )
    return s

if __name__ == '__main__':
    args = parse_cli()

    if args.nr:
        args.output = 'nr'
    threads = args.threads
    jobs = []
    n = len(args.accession)//threads if len(args.accession)//threads >0 else 1
    sub_acc = [args.accession[x:x+n] for x in range(0, len(args.accession), n)]

    for x in range(len(sub_acc)):
        p = Process(target = convert, args = (sub_acc[x], args.output))
        jobs.append(p)
        p.start()
    try:
        for p in jobs:
            p.join()
    except KeyboardInterrupt:
        for p in jobs:
            p.terminate()
    if args.fun:
        fun()
    #RefSeq

