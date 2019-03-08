#!/usr/bin/env python3

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.core.cli as corecli
import rotifer.core.functions as cf
import rotifer.neighborhood.neighborhood as nh
import rotifer.table.table as tb
import argparse
import warnings
import numpy as np
warnings.filterwarnings('ignore')
from collections import defaultdict
import subprocess
import argcomplete
import re
#TODO: Build class for input, output, filter

__version__ = 0.18
__authors__ = 'Gilberto Kaihami'

__doc__ = '''
NAME
rnexplorer

SYNOPSIS
Annotate, clean, reformat rneighbors output
Another line

DESCRIPTION

AUTHOR
'''

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Annotate, clean, reformat rneighbors output',
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
                        help = '''Add annotation column (usage: FILE/COL_KEY, COL_KEY:FILE, or COL_KEY:FILE:COL_NAME).''',
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

    parser.add_argument('-fr', '--filter_row',
                        dest = 'filter_row',
                        default = [],
                        action = 'append',
                        help = '''Custom filter operating over a row.
"COL_KEY:lambda row: expression" or "lambda row: expression".
USE QUOTES'''
                        )

    parser.add_argument('--doc',
                        help = 'Show documentation (--doc {[md/mardown], raw, man}',
                        dest = 'doc',
                        action = 'store_true')

    parser.add_argument('-oa', '--oa',
                        '--outputargs',
                        dest = 'outputargs',
                        default = [],
                        action = 'append',
                        help = 'Advanced option arguments'
                        )

    parser.add_argument('-v', '--verbose',
                        dest = 'verbose',
                        action = 'store_true',
                        help = 'Verbose')

    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rnexplorer',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Convert neighborhood output (table<->gi2operons)'
                                                  ))

    parser2 = corecli.config().input()

    parser_merged = corecli.merge_parser(parents = [parser, parser2], add_help = True)

    argcomplete.autocomplete(parser_merged)

    args = corecli.parseargs(parser_merged)

    return args

def exclude_type(df, exclude_dict):
    for k,v in exclude_dict.items():
        df = df[~(df[k].isin(v))]
    return df

def include_type(df, include_dict):
    for k,v in include_dict.items():
        df = df[df[k].isin(v)]
    return df

def sort_blocks_by(df, condition):
    #TODO finish this
    g = df.groupby('block_id').agg(dc)

def custom_row_filter(df, ipt_parameter, verbose = False,
                      new_info_ls = []):
    '''
    Custom function to operate over a row using a lambda function
    This function use pandas apply function for each row.
    Therefore the lambda expression should be like:
    lambda row: row[col] if ... else...
    For more complex operations check exec function
    -------------
    INPUT:
    df:            a dataframe
    ipt_parameter: a list containing the lambda expression
                   and column in the format (lambda_exp:column_name)
    verbose:       a boolean true or false
    new_info_ls:   a list containing columns
    -------------
    OUTPUT:
    A tuple containing a list of new columns and a df with the new added columns
    (new_info_ls, df)
    -------------
    USAGE:
    ipt_parameter = ['new_col1:lambda x: x**2',
                     'new_col2:lambda x: x/2']
    lambda row: row['col_name'] ...
    lambda row: row[col_number] ...

    '''
    ls_of_operations = []
    for e in ipt_parameter:
        _ = e.split(':') # check if is col_name:expression or expression only
        if len(_) <=2:
            col = len(df.columns)
            exp = ':'.join([x.strip() for x in _])
            if verbose:
                cf.vmsg(['Applying custom filter',
                         'Output column {0}'.format(col),
                         'Expression    "{0}"'.format(exp)])

            df[str(col)] = df.apply(eval(exp), axis = 1)
            new_info_ls.append(str(col))
        else:
            col = _[0]
            exp = ':'.join([x.strip() for x in _[1:]])
            if verbose:
                cf.vmsg(['Applying custom filter',
                         'Output column {0}'.format(col),
                         'Expression    "{0}"'.format(exp)])

            df[str(col)] = df.apply(eval(exp), axis = 1)
            new_info_ls.append(str(col))

    return (new_info_ls, df)

def exec_function(df, s, verbose = False):
    '''
    Execute a complex function.
    Class, functions object attributions can be done using this.
    '''
    try:
        func = open(s).read()
    except:
        func = s

    exec(s)
    return df


def strlist2dict(ls_of_str, output = 'str'):
    '''
    convert a list with strings in the format : or =
    The input could be:
    'KEY':'VALUE' or 'KEY'
    if only key is passed the value is set to 1
    ----
    select if output is str or list

    '''
    if output == 'str':
        dc = defaultdict(str)
        for element in ls_of_str:
            splitted = re.split('[:|=]', element)
            if len(splitted) == 1:
                dc[splitted[0]] = 1
            if len(splitted) == 2:
                dc[splitted[0]] = splitted[1]
    else:
        dc = defaultdict(list)
        for element in ls_of_str:
            splitted = re.split('[:|=]', element)
            if len(splitted) == 1:
                dc[splitted[0]].append(1)
            if len(splitted) == 2:
                dc[splitted[0]].append(splitted[1])

    return dc


def arg2dict(ls_of_str):
    dc = defaultdict(list)
    for select in ls_of_str:
        k,v = re.split('[:|=]', select)
        try:
            a = open(v).read().splitlines()
            dc[k].extend(a)
        except:
            dc[k].append(v)
    return dc

def ann2dict(ls_of_str, df):
    dc = defaultdict(list)
    for select in ls_of_str:
        splitted = re.split('[:|=]', select)
        if len(splitted) == 1: # map to pid and its a file or it is a column
            try:
                _ = pd.read_csv( splitted[0], header = None, sep = "\t")
                if len(_.columns) == 2:
                    _.columns = ['source', 'target']
                else:
                    s,t, *other = list(_.columns)
                    _.columns = ['source', 'target'] +other

                dc['pid'].append(_[['source', 'target']])


            except:
                try:
                    #include_col = tb.select_col(df.columns, args.includecol)
                    n = tb.select_col(df.columns, [splitted[0]])
                    for e in n:
                        if e in df.columns: #tb.select_col(df.columns, args.includecol)
                            dc['add_col'].append(e)
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
        if k != 'add_col':

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

if __name__ == '__main__':
    '''
    Variables declared
    +-----------------+--------------------+---------------------------------------------------------------+
    |      name       |    args            | why                                                           |
    +-----------------+--------------------+---------------------------------------------------------------+
    | args            | parser_cli()       | Parse user arguments                                          |
    | verbose         | args.verbose       | Check if verbose is activated                                 |
    | outputargs      | defaultdict()      | A dictionary contaning custom outputargs                      |
    | outformat       | args.outformat     | Check if user passed an output format otherwise same as input |
    | sym_path        |                    | Get the rotifer path                                          |
    | doc             |                    | Path to the documentation (man)                               |
    | s               | subprocess         | A process to open the doc                                     |
    | padtable        | args.padtable      | Check if pad is on                                            |
    | ipt             | nh.ipt2table       | Parse the input to a dataframe and get the input format       |
    | df              | ipt.df             | The input in a pandas dataframe format                        |
    | other_info      | args.addinfo       | Other info to put in the gi2operon header                     |
    | excludebytype   | args.excludebytype | Which rows should be excluded                                 |
    | includebytype   | args.includebytype | Which rows should be included                                 |
    | exclude_dict    | defaultdict()      | A dictionary containing the columns to exclude                |
    | include_dict    | defaultdict()      | A dictionary containing the columns to include                |
    | exclude_col     | tb.select_col()    | A list containing the columns to be excluded                  |
    | include_col     | tb.select_col()    | A list containing the columns to be included                  |
    | annotation      | args.annotation    | A raw input containing the annotations to be passed           |
    | annotation_dict | ann2dict()         | A dictionary containing the new information to be added       |
    | new_info_ls     | list               | A list containing all new information columns to be added     |
    +-----------------+--------------------+---------------------------------------------------------------+
    '''

    args = parse_cli()

    # verbose
    verbose = args.verbose

    ### parser optionargs
    outputargs = defaultdict(str)
    if args.outputargs:
        outputargs = strlist2dict(args.outputargs)

    outformat = args.outformat.lower()
    # Check if user called doc
    if args.doc:
        doc_format = ''.join(args.file)
        if  doc_format == 'man':

            sym_path = os.path.abspath(os.path.join(os.path.realpath(__file__), '../..'))

            doc = os.path.join(sym_path, 'doc/man/{0}'.format(os.path.basename(__file__).split('.')[0]))
            s = subprocess.Popen(['man {0}.man'.format(doc)], shell = True)
            s.communicate()[0]
            sys.exit()
        else:

            mdown = cf.rdoc(os.path.realpath(__file__))
            mdown.writer(doc_format)
            sys.exit()


    padtable = args.padtable

    if outformat:
        if outformat in ['gi2operon','gi2operons']:
            outformat = 'gi2operon'
        elif outformat in ['table', 'tables']:
            outformat = 'table'
        elif outformat in ['compact']:
            outformat = 'compact'
        elif outformat in ['gi2']:
            outformat = 'gi2'
        else:
            sys.stderr.write('Unknown option {0}\n'.format(args.outformat))
            sys.stderr.write('Select a valid output format')
            sys.exit()
    # ipt
    if verbose:
        cf.vmsg('Reading table')

    ipt = nh.ipt2table(args.file, verbose)
    df = ipt.df

    if verbose:
        cf.vmsg('Table loaded')

    # Custom

    # sort_blocks_by(df, ['EME40621.1', 'BAA32817.1'])

    # Additional information

    other_info = args.addinfo # Gi2operon header

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
        annotation_dict = ann2dict(annotation, df)

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

    new_info_ls = []
    if annotation_dict:
        new_info_ls, df = add_annotation(df, annotation_dict)

    exclude_col = []

    if args.excludecol:
        exclude_col = tb.select_col(df.columns, args.excludecol)

    include_col = []
    if args.includecol:
        include_col = tb.select_col(df.columns, args.includecol)

    if args.filter_row:
        new_info_ls, df = custom_row_filter(df, args.filter_row,
                                            new_info_ls = new_info_ls,
                                            verbose = verbose)

    # Add additional columns
    if 'add_col' in annotation_dict.keys():
        for ele in annotation_dict['add_col']:
            new_info_ls.append(ele)

    if verbose:
        cf.vmsg('Starting to print')

    if not df.empty:
        nh.writer(df, of = outformat, other_info = other_info,
                   new_info_ls=new_info_ls, pad = padtable,
                   exclude_col=exclude_col,
                  include_col=include_col, verbose = args.verbose, kwargs = outputargs)
