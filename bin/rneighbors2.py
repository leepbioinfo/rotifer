#!/usr/bin/env python3

__version__ = '0.51'
import os
import sys

# Rotifer modules
import rotifer.core.cli as corecli
from rotifer.db.neighbors import rneighbors
from rotifer.db.cli import rneighbors as clineighbors
import rotifer.db.cli as sqlcli
# Other imports
import argparse # Not our command line parser of choice but better than none

def parse_cli():

    parse = sqlcli.rneighbors().input()
    parse2 = sqlcli.sql().input()
    args = corecli.parseargs(parents = [parse, parse2])
    return args

if __name__ == '__main__':
    args = parse_cli()
    neighbors.rneighbors(above = (args.above),
                         below = (args.below),
                         outformat = args.outformat,
                         input_file = (args.file),
                         gacc = args.genomic_accession,
                         asm = args.assembly,
                         database = (args.database),
                         port = (args.port),
                         user = args.user,
                         debug = args.debug,
                         log = args.log
                         ).runRneighbors()
