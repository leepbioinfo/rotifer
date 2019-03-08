#!/usr/bin/env python3
import os
import sys
sys.path.insert(0,'/home/kaihami/mymodules')
import rotifer.core.cli as corecli
import argparse
def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'accession',
        help = "Input accession file",
        nargs = '*',
        action = corecli.action.autoload, duplicates = False
        )
    args = parser.parse_args()
    return args

args = parse_cli()

data = [x.split('\t') for x in args.accession]
col_width = max(len(word) for row in data for word in row) + 2  # padding
for row in data:
    print ("".join(word.ljust(col_width) for word in row))
