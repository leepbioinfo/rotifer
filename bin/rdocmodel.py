#!/usr/bin/env python3

import sys
import os
import rotifer.core.functions as cf
import subprocess
import argparse
import mdv
__version__ = 0
__authors__ = ''

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Build a model for docs (markdown format)')
    parser.add_argument('file',
                        help = 'Insert file'
                        )
    parser.add_argument('-o',
                        '--output',
                        help = 'Output format (markdown, man, raw)',
                        default = 'markdown')

    args = parser.parse_args()

    return args
def fun(f):
    toprint = toprint.rstrip('\n')

    if mkdocs_output == 'raw':
        print(toprint)
        sys.exit()
    if mkdocs_output in ['md','markdown']:
        formatted = mdv.main(toprint, theme='641.4291')
        print(formatted)

    if mkdocs_output in ['man']:
        pass


if __name__ == '__main__':
    args2 = parse_cli()
    mkdocs_output = args2.output
    mark = cf.rdoc(args2.file)
    mark.writer(mkdocs_output)

