#!/usr/bin/env python3
import os
import sys
import argparse
sys.path.insert(0, '/home/kaihami/mymodules')

import rotifer.core.cli as corecli

import argcomplete
from tempfile import mkstemp
import subprocess

__version__ = '0.2'
__authors__ = 'Gilberto Kaihami; Aureliano Guedes'
__rdoc__='''
DESCRIPTION:
Get a distance between PDBs
'''


def parse_cli():
    parser = corecli.parser(description = 'Get a distance between PDBs')

    parser.add(dest = 'list',
                        helper = 'List of PDBs',
                        nargs = '*',
                        action = corecli.action.autoload)
    parser.add(short_arg = '-a',
               long_arg = '--align',
                        helper = 'Alignment types (Default: align, super, cealign)',
                        action = 'append',
                        default = [])
    parser.add(short_arg = '-t',
               long_arg = '--threads',
                        helper = 'Number of threads (default = 3)',
                        arg_type = str,
                        default = '3')

    argcomplete.autocomplete(parser)
    args = parser.parse_args() 

    return args

if __name__ == '__main__':
    args = parse_cli()
    fd, path_pdbs = mkstemp()
    if not args.align:
        align = ['align', 'super', 'cealign']
    else:
        align = [x.lower() for x in args.align]
    with open(path_pdbs, 'a') as f:
        f.write('\n'.join(args.list))
    try:
        p = subprocess.Popen(['pymol -cQ /home/kaihami/mymodules/rotifer/bin/pymol_rmsd_comparer.py {0} {1} '.format(path_pdbs, args.threads) +' '.join(align)], shell = True)
        p.communicate()
    except KeyboardInterrupt:
        sys.exit(0)
