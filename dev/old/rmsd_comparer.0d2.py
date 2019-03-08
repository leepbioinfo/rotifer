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


def parse_cli():
    parser = argparse.ArgumentParser(description = 'Get a distance between PDBs',
                                     formatter_class = argparse.RawTextHelpFormatter)

    parser.add_argument('list',
                        help = 'List of PDBs',
                        nargs = '*',
                        action = corecli.action.autoload)
    parser.add_argument('-a', '--align',
                        help = 'Alignment types (Default: align, super, cealign)',
                        action = 'append',
                        default = [])
    parser.add_argument('-t', '--threads',
                        help = 'Number of threads (default = 3)',
                        type = str,
                        default = '3')
    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rmsd_comparer',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Input a list of PDBs return the RMSD between all PDBs'))
    parser2 = corecli.config().input()

    parser_merged = corecli.merge_parser(parents = [parser, parser2],
                                         add_help = False)
    argcomplete.autocomplete(parser_merged)
    args = corecli.parseargs(parser_merged)

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
