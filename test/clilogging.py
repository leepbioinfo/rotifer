#!/usr/bin/env python3

import sys
import rotifer
import rotifer.core.cli as corecli

__version__ = 0.001
__authors__ = 'YOUR NAME'

def parse_cli():
    parser = corecli.parser(description = 'PROGRAM DESCRIPTION')
    parser.add( long_arg = '--verbose',
                short_arg = '-v',
                dest = 'verbose',
                helper = 'Some help',
                action = "count"
                )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_cli()
    verbose = args.verbose
    if verbose:
        rotifer.logger.info(f'{sys.argv}')
        rotifer.logger.warning('This is warning')
        rotifer.logger.debug('This is debug')
