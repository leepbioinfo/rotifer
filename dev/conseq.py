#!/usr/bin/env python3

import os
import argparse
import sys
sys.path.insert(0, '/home/kaihami/mymodules')

import rotifer.core.cli as corecli
import argcomplete
import pandas as pd

__version__ = 0.3
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser(description = 'Extract a column from a sequence')
    parser.add(
               dest = 'fasta')

    parser.add(short_arg= '-p',
               long_arg = '--position',
               dest = 'position',
               action='append',
               helper = 'Position in the reference sequence (number in the reference sequence or a string)')

    parser.add(short_arg = '-c',
               long_arg = '--count',
               dest = 'count',
                        helper = 'Count frequency based on the reference accession',
                        action = 'store_true')

    parser.add(short_arg = '-a',
               long_arg = '--accession',
               dest = 'accession',
               default = '',
               helper = 'Reference Sequence')

    parser.add(short_arg = '-of',
               long_arg = '--outformat',
               dest = 'format',
               helper = 'Output format [String, Table]',
               default = 'table')

    parser.add(long_arg = '--postion_type',
               helper = 'Select if position is a number or string',
               default = 'int')

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    return args
def fasta2df(fi):
    tmp = {}
    for x in range(0, len(fi)):
        if fi[x].startswith('>'):
            header = fi[x]
            header = fi[x].replace('>','')
            tmp[header] = ''
        else:
            tmp[header] += fi[x]
    df = pd.DataFrame({'Header': [x.split(' ')[0] for x in tmp.keys()],
                       'Seq': list(tmp.values())})
    return df

def count_matches(x, ls):
    tot = 0
    for f,s in zip(x,ls):
        if f == s:
            tot +=1
    return str(tot)

if __name__ == '__main__':
    args = parse_cli()
    df = fasta2df(open(args.fasta).read().splitlines())
    s = {}
    refseq = ''.join([x for x in df[df['Header'].str.contains(args.accession, regex = True)].Seq.values[0] if x != '-'])
    list_pos = []
    original = df[df['Header'].str.contains(args.accession, regex = True)].Seq.values[0]
    for pos in args.position:
        if '..' in pos:
            for pos2 in range(int(pos.split('..')[0]),
                              int(pos.split('..')[1])+1) :
                pos_refseq = 0
                pos_original = 0
                for e in original:
                    if e != '-':
                        if pos_original == int(pos2)-1:
                            if original[pos_refseq]+str(pos_original+1) not in s.keys():
                                s[original[pos_refseq]+str(pos_original+1)] = pos_refseq
                                list_pos.append(original[pos_refseq])
                        pos_original +=1
                    pos_refseq +=1
        else:
            pos_refseq = 0
            pos_original = 0
            for e in original:
                if e != '-':
                    if pos_original == int(pos)-1:
                        if original[pos_refseq]+str(pos_original+1) not in s.keys():
                            s[original[pos_refseq]+str(pos_original+1)] = pos_refseq
                            list_pos.append(original[pos_refseq])
                    pos_original +=1
                pos_refseq +=1

    for k,v in s.items():
        df[k] =df['Seq'].map(lambda x: x[v])

    for k in s.keys():
        if 'sub' in df.columns:
            df['sub'] = df['sub'] +df[k]
        else:
            df['sub'] = df[k]
    df['count'] = df['sub'].map(lambda x: count_matches(x, list_pos))
    if args.format.lower() == 'table':
        if args.count:
            print('\t'.join(['Accession'] + [k for k in s.keys()] + ['count']))
        else:
            print('\t'.join(['Accession'] + [k for k in s.keys()]))
        for i,row in df.iterrows():
            if args.count:
                cols = ['Header'] + list(s.keys()) +['count']
                print('\t'.join([row[col] for col in cols]))
            else:
                cols = ['Header'] + [k for k in s.keys()]
                print('\t'.join([row[col] for col in cols]))
    else:
        _ = ''.join([k for k in s.keys()])
        if args.count:
            print('\t'.join(['Accession',_, 'count'])) 
        else:
            print('\t'.join(['Accession',_])) 
        for i,row in df.iterrows():
            seq2write = ''.join([row[k] for k in s.keys()])
            if args.count:
                print('\t'.join([row['Header'], seq2write, row['count']]))

            else:
                print('\t'.join([row['Header'], seq2write]))



