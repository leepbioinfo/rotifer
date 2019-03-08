#!/usr/bin/env python3

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.core.cli as corecli
import argparse
from multiprocessing import Process
import warnings
import numpy as np
warnings.filterwarnings('ignore')

__version__ = 0.02
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Convert rneighbors output (table<->gi2operons)',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('file',
                        help = 'Input neighborhood file',
                        action = corecli.action.autoload,
                        nargs = '*',
                        duplicates = False)

    parser.add_argument('-y', '--header',
                              action='append',
                              help = 'Add more information to gi2operon header. The additional information must be a valid collumn in the table format\n\
(example: -y assembly, -y classification -y seq_type)',
                              default = [])

    parser.add_argument('-of', '--outformat',
                        help = 'Output format (table, gi2ioperon, not implemented)',
                        default = '')

    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rnexplorer',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Convert neighborhood output (table<->gi2operons)'
                                                  ))

    args = parser.parse_args()

    return args

def collen(df):
    df['start'].astype(str)
    df['end'].astype(str)
    df['plen'].astype(str)
    q = '-->'
    df['cds_loc'] = df['start'].astype(str)+'..'+df['end'].astype(str)
    direction = len('dir')
    len_plen = df['plen'].astype(str).str.len().max()
    len_pid = df['pid'].str.len().max() if df['pid'].str.len().max() >= len('pid') else len('pid')
    len_type = df['type'].str.len().max() if df['type'].str.len().max() >= len('type') else len('type')
    len_gene = df['gene'].str.len().max() if df['gene'].str.len().max() >= len('gene') else len('gene')
    len_cds = df['cds_loc'].str.len().max()
    len_locus = df['locus'].str.len().max() if df['locus'].str.len().max() >= len('locus') else len('locus')
    try:
        len_modified = df['modified'].str.len().max() if df['modified'].str.len().max() >= len('gi') else len('gi')
    except:
        len_modified = 1
    len_product = 0
    collen_len = [len(q), len_cds, len('dir'), len_plen, len_pid, len_type,
                  len_gene, len_locus, len_modified, len_product]

    return collen_len

def gi2operon(sub_df, other_info = []):

    sub_df['plen'].astype(str)
    if '' in sub_df['plen'].values:
        sub_df.loc[sub_df[sub_df['plen'] == ''].index, 'plen'] = '.'

    if '' in sub_df['pid']:
        sub_df.loc[sub_df[sub_df['pid'] == ''].index, 'pid'] = '.'

    col_len = collen(sub_df)
    protein = sub_df[sub_df['query'] == '1']['pid'].values[0]

    nucleotide = sub_df.nucleotide.values[0]

    max_header = 'ORGANISM ' + df.organism.values[0] + ' accession no is ' + nucleotide + ' Protein is '+ protein
    if other_info:
        for e in other_info:
            max_header += ' | '+ e + ':'+''.join(str(sub_df[e].values[0]))

    print(max_header)

    header = [".","cds","dir","len","pid","type","gene","locus","gi",'product']
    header = [header[int(i)].ljust(int(col_len[int(i)])) for i in range(0,len(header), 1)]
    print('  '.join(header))

    # sub_df['.'] = np.where(sub_df['query'] == 1, '-->', '.')
    # sub_df['dir'] = np.where(sub_df['strand'] == 1, '+', '-')
    # sub_df['cds'] = sub_df['start'].astype(str) + '..' + sub_df['end'].astype(str)
    # sub_df['len'] =sub_df['plen']
    # sub_df['gi'] = sub_df.shape[0]*['.']
    # print(sub_df.to_string(columns = header, index = False, justify = 'left'))
    #
    for i, row in sub_df.iterrows():
        query = '-->' if row.pid in sub_df[sub_df['query'] == 1].pid.values else '.'
        direction = '+' if row.strand == '1' else '-'

        cds = str(row.start) +'..'+str(row.end)
        toprint = [str(x) for x in [query, cds, direction,
                          row.plen, row.pid, row.type,
                          row.gene, row.locus, '.', row['product']]]
        toprint = [toprint[i].ljust(int(col_len[i])) for i in range(len(toprint))]
        print('  '.join(toprint))
    print("---------------------------------------")

def table(f):
    ipt = f

    block_id = 0
    header = 'nucleotide start end strand block_id, query pid type plen locus seq_type assembly gene product organism classification'.split()

    print('\t'.join(header))

    for block in ipt:
        if block != '':
            block_id +=1
            block_split = [x for x in block.split('\n') if x != ''and x != "---------------------------------------"]

            for x in range(len(block_split)):
                if x == 0:
                    chr_type = '.'
                    asm = '.'
                    classification = '.'

                    organism = block_split[x].split(' accession no')[0].strip()
                    nucleotide = block_split[x].split(' accession no is')[1].split('Protein')[0].strip()

                    if '|' in block_split[x]:
                        other_information = block_split[x].split('|')

                        for ele in other_information:
                            if 'seq_type' in ele:
                                chr_type = ele.split(':')[1].strip()
                            if 'assembly' in ele:
                                asm = ele.split(':')[1].strip()
                            if 'classification' in ele:
                                classification = ele.split(':')[1].strip()

                if x>1:
                    line = [x for x in block_split[x].split(' ') if x != '']
                    query, start_end, strand, length, pid,cds_type, gene, locus, gi, *product = line
                    start, end = start_end.split('..')
                    query = 1 if query == '-->' else 0
                    strand = 1 if strand == '+' else -1
                    print('\t'.join([nucleotide, start,end, str(strand), str(block_id), str(query),
                                     pid, cds_type, length, locus, chr_type, asm, gene,
                                     ' '.join(product), organism, classification]))


if __name__ == '__main__':
    args = parse_cli()
    # to_open = open(args.file)
    to_open = args.file

    other_info = args.header

    if args.outformat:
        pass

    else:
        if 'ORGANISM' in to_open[0]:
            format_type = 'gi2operon'

        else:
            format_type = 'table'

        if format_type == 'table':
            # print(to_open)
            s = [x.split('\t') for x in to_open]

            df = pd.DataFrame(s[1:], columns = s[0])

            df = df.fillna('.')
            g = df.groupby('block_id')
            for i, group in g:
                group = group.fillna('.')
                gi2operon(group, other_info)

        if format_type == 'gi2operon':
            table(to_open)

