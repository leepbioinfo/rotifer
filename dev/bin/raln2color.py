#!/usr/bin/env python3

from rotifer.table import table as tb

import rotifer.core.cli as corecli

__version__ = 0.001
__authors__ = 'Gilberto Kaihami; Robson F. Souza'

# TODO: add consensus
# TODO: add other colors
# TODO: show residues by above some threshold

def parse_cli():
    parser = corecli.parser(description = 'Add color to alignment')

    # Add another options here

    parser.add(
                dest = 'fasta',
                nargs = '*',
                helper = 'Fasta file',
                action = corecli.action.autoload,
                duplicates = False
                )

    parser.add( long_arg = '--scale',
                short_arg = '-s',
                dest = 'scale',
                helper = 'Add a scale row',
                action = "store_true"
                )

    parser.add( long_arg = '--colorscheme',
                short_arg = '-cs',
                dest = 'colorscheme',
                default = 'Clustal',
                arg_type = str,
                helper = 'Select color scheme (Clustal)',
                action = "store"
                )

    parser.add( long_arg = '--color',
                short_arg = '-c',
                dest = 'color',
                default = 'background',
                arg_type = str,
               helper = 'Select color (default: background) [background/foreground]',
                action = "store"
                )

    args = parser.parse_args()

    return args

def color_bg(s, color = ''):
    '''
    s: String
    '''
    if color:
        color = f'48;5;{color}'

    return f'\033[{color}m{s}\033[m'

def color_fg(s, color = ''):
    if color:
        color = f'38;5;{color}'

    return f'\033[{color}m{s}\033[m'

def color_res(s, cs):
    if s in 'A I L M F W V'.split():
        return cs(s, 33)

    if s in 'K R'.split():
        return cs(s, 124)

    if s in 'E D'.split():
        return cs(s, 127)

    if s in 'N Q S T'.split():
        return cs(s, 34)

    if s  == 'C':
        return cs(s, 168)

    if s == 'G':
        return cs(s, 166)

    if s == 'P':
        return cs(s, 178)

    if s in 'H Y'.split():
        return cs(s, 37)

    else:
        return s

if __name__ == '__main__':
    args = parse_cli()

    df = tb.fasta2df(args.fasta) # seq is a list

    select = {'background':color_bg,
              'bg':color_bg,
              'foreground': color_fg,
              'fg': color_fg}

    if args.scale:
        scale = len(df['Seq'].values[0])
        scale_dot = ''
        scale_number = ''
        scale_bar = ''
        for x in range(0, scale):
            if x == 0:
                scale_number = str(1)
                scale_bar += '|'
                scale_dot += '-'
            if x % 10 == 0:
                scale_number =scale_number[:-1]

                scale_number += str(x)
                scale_number += ' '*(10-len(str(x+1)))
                scale_number += ' '
                scale_bar = scale_bar [:-1]

                scale_bar += '|'
                scale_bar += ' '
                scale_dot += '-'
            else:
                scale_dot += '-'
                scale_bar += ' '

    df['colored'] = df['Seq'].map(lambda x: ''.join([color_res(y, select[args.color]) for y in x]))
    # col_len

    col_len = max([len(x) for x in df['ID'].values])

    df['ID_len'] = df['ID'].map(lambda x: x.ljust(col_len))
    df['2print'] = df['ID_len'] + '  ' + df['colored']


    if args.scale:
        print(' '* col_len + '  '+ scale_number)
        print(' '* col_len + '  '+ scale_bar)
        print(' '* col_len + '  '+ scale_dot)
    print('\n'.join(df['2print'].values))

