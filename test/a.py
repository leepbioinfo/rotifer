#!/usr/bin/env python3

import rotifer.core.cli as corecli
from rotifer.core.log import log
__version__ = 0.001
__authors__ = 'YOUR NAME'

def parse_cli():
    parser = corecli.parser(description = 'PROGRAM DESCRIPTION')

    # Add another options here

    parser.add( long_arg = '--verbose',
                short_arg = '-v',
                dest = 'verbose',
                helper = 'Some help',
                action = "count"
                )

    parser.add( long_arg = '--log_file',
                dest = 'log_file',
                default = '',
                arg_type = str,
                helper = 'Some help',
                action = "store"
                )

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()
    verbose = args.verbose

    if verbose:
        import sys
        log({1:f'{sys.argv}',
             2:'This is warning',
             3:'This is debug'}, level = verbose,
            log_file = args.log_file)
