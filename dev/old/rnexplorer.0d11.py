#!/usr/bin/env python3

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.core.cli as corecli
import rotifer.neighborhood.neighborhood as nh
import rotifer.table.table as tb
import argparse
import warnings
import numpy as np
warnings.filterwarnings('ignore')
from collections import defaultdict
import subprocess

#TODO: Build class for input, output, filter

__version__ = 0.11
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Annotate, clean, reformat rneighborsoutput',
                                     formatter_class = argparse.RawTextHelpFormatter,
                                     add_help = False)
    parser.add_argument('file',
                        help = 'Input neighborhood file',
                        # required='--configdump' not in sys.argv,
                        action = corecli.action.autoload,
                        nargs = '*',
                        duplicates = False)

    parser.add_argument('-y', '--header',
                        action='append',
                        dest = 'addinfo',
                              help = 'Add more information to gi2operon header. The additional information must be a valid collumn in the table format\n\
(example: -y assembly, -y classification -y seq_type)',
                              default = [])

    parser.add_argument('-of', '--outformat',
                        dest = 'outformat',
                        help = 'Output format (table, gi2ioperon)',
                        default = '')

    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rnexplorer',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Convert neighborhood output (table<->gi2operons)'
                                                  ))
    parser.add_argument('--xt', '--exclude_by_type',
                        '-exclude_by_type', '-xt',
                        dest = 'excludebytype',
                        default = [],
                        help = '''Exclude individual features type (usage: KEY:VALUE).''',
                        action = 'append',
                        type = str)

    parser.add_argument('--it','--include_by_type',
                        '-include_by_type', '-it',
                        dest = 'includebytype',
                        default = [],
                        help = '''Include features type (usage: KEY:VALUE).''',
                        action = 'append',
                        type = str)
    parser.add_argument('-a', '--annotation',
                        '--a', '-annotation',
                        dest = 'annotation',
                        default = [],
                        help = '''Add annotation column (usage: FILE, COL_KEY:FILE, or COL_KEY:FILE:COL_NAME).''',
                        action = 'append',
                        type = str)

    parser.add_argument('-p', '--padtable',
                        dest = 'padtable',
                        help = 'Format gi2operon output',
                        action = 'store_true')

    parser.add_argument('-xc', '--exclude_column',
                        '--xc', '-exclude_column',
                        dest = 'excludecol',
                        help = '''Exclude a column from output (usage: col_number, start..end, column name).''',
                        action = 'append',
                        type = str,
                        default = [])

    parser.add_argument('-ic', '--include_column',
                        '--ic', '-include_column',
                        dest = 'includecol',
                        help = '''Include a column from output (usage: col_number, start..end, column name).''',
                        action = 'append',
                        type = str,
                        default = [])
    parser.add_argument('-ai','--add_information',
                        '--ai',
                        help = '''Include an additional column in the gi2operon.
example: -ai arch''',
                        dest = 'additional',
                        action = 'append',
                        type = str,
                        default = [])

    parser.add_argument('-c',
                        '--compact',
                        dest = 'compact',
                        default = '',
                        type = str,
                        help = 'Output compact format (usage: -c architecture_column')

    parser.add_argument('--doc',
                        help = 'Show documentation (--doc 0)',
                        dest = 'doc',
                        action = 'store_true')

    parser2 = corecli.config().input()
    args = corecli.parseargs(parents = [parser, parser2], exclude_from_dump = ['file'])

    return args

def exclude_type(df, exclude_dict):
    for k,v in exclude_dict.items():
        df = df[~(df[k].isin(v))]
    return df

def include_type(df, include_dict):
    for k,v in include_dict.items():
        df = df[df[k].isin(v)]
    return df

def arg2dict(ls_of_str):
    dc = defaultdict(list)
    for select in ls_of_str:
        k,v = select.split(':')
        try:
            a = open(v).read().splitlines()
            dc[k].extend(a)
        except:
            dc[k].append(v)
    return dc

def ann2dict(ls_of_str):
    dc = defaultdict(list)
    for select in ls_of_str:
        splitted = select.split(':')
        if len(splitted) == 1: # map to pid and its a file
            try:
                # a = open(splitted[0]).read().splitlines()
                # source_ls = []
                # target_ls = []
                # for line in a:
                #     s = line.split('\t')
                #     source_ls.append(s[0])
                #     target_ls.append(s[1])
                _ = pd.read_csv( splitted[0], header = None, sep = "\t")
                if len(_.columns) == 2:
                    _.columns = ['source', 'target']
                else:
                    s,t, *other = list(_.columns)
                    _.columns = ['source', 'target'] +other

                dc['pid'].append(_[['source', 'target']])


            except:
                sys.stderr.write('Unable to read {0}\n'.format(splitted[0]))
                sys.exit()

        if len(splitted) == 2: # could be pid:file  or column:file
            try:
                _ = pd.read_csv( splitted[1], header = None, sep = "\t")
                if len(_.columns) == 2:
                    _.columns = ['source', 'target']
                else:
                    s,t, *other = list(_.columns)
                    _.columns = ['source', 'target'] +other

                dc[splitted[0]].append(_[['source', 'target']])
            except:
                sys.stderr.write('Unable to read {0}\n'.format(splitted[1]))
                sys.exit()
        if len(splitted) == 3: # collumn:file:column_name
            try:
                _ = pd.read_csv( splitted[1], header = None, sep = "\t")
                if len(_.columns) == 2:
                    _.columns = ['source', 'target']
                else:
                    s,t, *other = list(_.columns)
                    _.columns = ['source', 'target'] +other

                dc[splitted[0]].append((_[['source', 'target']], splitted[2]))
            except:
                sys.stderr.write('Unable to read {0}\n'.format(splitted[1]))
                sys.exit()

    return dc


def add_annotation(df, dc):
    new_info_ls = []
    for k in dc.keys(): # k is the column

        for ann_df in dc[k]: # each df
            if len(ann_df) !=2:
                try: # try to merge
                    df = df.merge(ann_df, how = 'left', left_on = k, right_on = 'source')
                    df.drop(columns= ['source'], inplace = True)
                    df.columns = (list(df.columns[:-1])+[str(len(df.columns))])
                    df = df.fillna('.')
                    new_info_ls.append(str(len(df.columns)))
                except: pass

            if len(ann_df) == 2:
                try: # try to merge
                    ann, col_ann = ann_df

                    df = df.merge(ann, how = 'left', left_on = k, right_on = 'source')
                    df.drop(columns= ['source'], inplace = True)
                    df.columns = (list(df.columns[:-1])+[col_ann])
                    df = df.fillna('.')
                    new_info_ls.append(col_ann)
                except: pass

    return (new_info_ls,df)
#     df[str(len(df.collumns))] #left join

def get_block_id(series):
    query = series.unique()
    if len(query) ==1:
        return query
    else:
        query = '; '.join(list(series.unique()))
        return query

def compact(df, col_arch= 'pfam', pad = True, new_info_ls = []):
    df = df.reset_index(drop = True)
    df = df.astype(str)
    df["block_id"] = pd.to_numeric(df["block_id"])
    if col_arch:
        df['arch_compact'] = np.where(df[col_arch] == '.', '?', df[col_arch])
        df['arch_compact'] = np.where(df['query'] == '1', df.arch_compact+'*',
                                       df.arch_compact)
        df['arch_compact'] = np.where(df.strand == '1', df['arch_compact']+'->', '<-' + df['arch_compact'])
        df['arch_compact'] = np.where( df['block_id'] == df.shift()['block_id'],
                                  np.where((df['strand'] == '1') & (df.shift().strand == '-1'),
                                                                     '||' +df['arch_compact'],
                                           df['arch_compact']),
                                           df['arch_compact'])
    dc = {'arch_compact':'sum'}
    if new_info_ls:
        for k in new_info_ls:
            dc[k] = get_block_id
    g = df.groupby('block_id').agg(dc).reset_index()
    if pad:
        g['block_id'] = g['block_id'].astype(str)
        g['left'] = g['arch_compact'].str.split('*').str[0]

        g['right'] = g['arch_compact'].map(lambda x: '*'+''.join(x.split('*')[1:] if len(x.split('*')) >1 else '*'   ))
        c_len = g['left'].str.len().max()
        cols2print = [x for x in new_info_ls if x not in ['block_id', 'right',
                                                   'left', 'arch_compact']]
        if cols2print:
            print('\t'.join(['block_id', 'neighborhood'] + cols2print))
        else:
            print('\t'.join(['block_id', 'neighborhood'] + cols2print))
        for i, row in g.iterrows():
            if cols2print:
                print('\t'.join([row.block_id]+
                                ['{0:>{1}}{2}'.format(row.left, c_len, row.right)]+
                                 list(row[cols2print].values)))
            else:
                print('\t'.join([row.block_id]+
                                ['{0:>{1}}{2}'.format(row.left, c_len, row.right)]))

    else:
        g['block_id'] = g['block_id'].astype(str)
        g['left'] = g['arch_compact'].str.split('*').str[0]

        g['right'] = g['arch_compact'].map(lambda x: '*'+''.join(x.split('*')[1:] if len(x.split('*')) >1 else '*'   ))
        cols2print = [x for x in new_info_ls if x not in ['block_id', 'right',
                                                   'left', 'arch_compact']]

        if cols2print:
            print('\t'.join(['block_id', 'neighborhood'] + cols2print))
        else:
            print('\t'.join(['block_id', 'neighborhood'] + cols2print))
        for i, row in g.iterrows():
            if cols2print:
                print('\t'.join([row.block_id]+
                                ['{0}{1}'.format(row.left, row.right)]+
                                 list(row[cols2print].values)))
            else:
                print('\t'.join([row.block_id]+
                                ['{0}{1}'.format(row.left, row.right)]))
if __name__ == '__main__':
    args = parse_cli()
    # print(args)
    outformat = args.outformat
    # sys.exit()
    # Check if user called doc
    if args.doc:
        sym_path = os.path.abspath(os.path.join(os.path.realpath(__file__), '../..'))

        doc = os.path.join(sym_path, 'doc/man/{0}'.format(os.path.basename(__file__).split('.')[0]))
        s = subprocess.Popen(['man {0}.man'.format(doc)], shell = True)
        s.communicate()[0]
        sys.exit()


    padtable = args.padtable

    if outformat:
        ofl = args.outformat.lower()
        if ofl.lower() in ['gi2operon','gi2operons']:
            outformat = 'gi2operon'
        elif ofl.lower() in ['table', 'tables']:
            outformat = 'table'
        else:
            sys.stderr.write('Unknown option {0}\n'.format(args.outformat))
            sys.stderr.write('Select a valid output format')
            sys.exit()
    # ipt
    ipt = nh.ipt2table(args.file)
    df = ipt.df
    new_info_ls = []

    # Additional information

    other_info = args.addinfo

    # Include/Exclude Block

    exclude_dict = defaultdict(list)
    include_dict = defaultdict(list)

    excludebytype = args.excludebytype
    if excludebytype:
        exclude_dict = arg2dict(excludebytype)

    includebytype = args.includebytype
    if includebytype:
        include_dict = arg2dict(includebytype)

    # Annotation Block

    annotation_dict = defaultdict(list)
    annotation = args.annotation

    if annotation:
        annotation_dict = ann2dict(annotation)

    # Select output format
    if outformat: pass
    else:
        outformat = 'gi2operon' if ipt.input_format_type == 'gi2operon' else 'table'

    # Apply exclude filter
    if exclude_dict and not include_dict:
        df = exclude_type(df, exclude_dict)

    # Apply include filter
    if include_dict:
         df = include_type(df, include_dict)

    # Apply annotation dict
    if annotation_dict:
        new_info_ls, df = add_annotation(df, annotation_dict)
    exclude_col = []

    if args.excludecol:
        exclude_col = tb.select_col(df.columns, args.excludecol)

    include_col = []
    if args.includecol:
        include_col = tb.select_col(df.columns, args.includecol)

    if args.additional:
        for ele in args.additional:
            new_info_ls.append(ele)

    if args.compact: # TODO
        compact(df, args.compact,
                pad = padtable,
                new_info_ls = new_info_ls)
#    print(include_col)


    #
    # if not df.empty:
    #     nh.printer(df, outformat, other_info,
    #                new_info_ls, padtable,
    #                exclude_col, include_col)
