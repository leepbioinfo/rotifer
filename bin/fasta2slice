#!/usr/bin/env python3

import os
import pandas as pd
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))

### Import core cli
import rotifer.core.cli as corecli
import argparse
import time
def parse_cli():
    parser = argparse.ArgumentParser(description='Extract sequence fragments from a fasta file')
    parser.add_argument('fasta',
                        help = 'Fasta file',
                        nargs = '*',
                        action = corecli.action.autoload, duplicates = False)
    parser.add_argument('-f', '--feature',
                        help = 'A feature table containing: Accession start end')
    args = parser.parse_args()
    return args

def fasta2df(fi):
    # Parse the input fasta file
    # Return: a dict with key = sequence header, value = list of nucleotides
    tmp = {}
    for x in range(0, len(fi)):
        if fi[x].startswith('>'):
            header = fi[x].replace('>','')
            tmp[header] = ''
        else:
            tmp[header] += fi[x]
    df = pd.DataFrame({'Header': [x.split(' ')[0] for x in tmp.keys()],
                       'rest': [' '.join(x.split(' ')[1:]) for x in tmp.keys()],
                       'Seq': list(tmp.values())})
    return df


if __name__ == '__main__':
    args = parse_cli()
    fasta_df = fasta2df(args.fasta)
    ft_df = pd.read_csv(args.feature, sep = '\t', header = None)
    ft_df.columns = ['acc', 'start', 'end']
    ft_df = ft_df.merge(fasta_df, left_on = 'acc', right_on = 'Header')
    ft_df.drop(columns = ['Header'], inplace = True)
    ft_df['number'] = ft_df.groupby('acc').cumcount()+1

    ft_df['partial_seq'] = ft_df.apply(lambda x: x['Seq'][x['start']-1:x['end']], 1)
    #print(ft_df.to_csv(sep = '\t'))
    ft_df['mod_header'] = ft_df['acc']+':'+ft_df['start'].astype(str)+'..'+ft_df['end'].astype(str)+':'+ ft_df['number'].astype(str)
    for i, v in ft_df.iterrows():
        print('>'+v.mod_header)
        print(v.partial_seq)
