#!/usr/bin/env python3

import rotifer.core.cli as corecli
from rotifer.table import table as tb
import pandas as pd

__version__ = 0.001
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser(description = 'Prepare data to ccmpred')

    # Add another options here
    parser.add(
                dest = 'fasta',
                nargs = '*',
                helper = 'Fasta file',
                action = corecli.action.autoload,
                duplicates = False
                )

    parser.add(
             '--reference',
             '-r',
             dest = 'reference',
             helper = 'Reference sequence id',
             arg_type = str)

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()

    df = tb.fasta2df(args.fasta)
    ref = args.reference

    # df = tb.fasta2df(open('amt.nr.0d8.amt.domain.sliced.filtered.formated.ordered.modified.trimmed.msa.iqtree.Rh.extended.msa').read().splitlines())
    # ref = '3HD6_A'

    n = df.Seq.map(lambda x: list(x))
    cols = n.to_frame().Seq.apply(pd.Series)

    list_pos = [x for x,y in enumerate(df[df['ID'] == ref].Seq.values[0]) if y != '-']
    cols = cols[list_pos].sum(1)
    cols = cols.to_frame()
    cols.columns = ['trimmed']
    df2 = df.join(cols)

    df2['2write'] = '>'+df2['Header'] +'\n'+df2['trimmed']
    print('\n'.join(df2['2write'].values))
