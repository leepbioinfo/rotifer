#!/usr/bin/env python3

import pandas as pd
from pandas.compat import StringIO
import warnings
import numpy as np
import sys
import os
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.core.functions as cf
from collections import defaultdict
import time

warnings.filterwarnings('ignore')
__version__ = 0.2
__authors__ = 'Gilberto Kaihami; Robson Souza'

class ipt2table:
    '''
    Input string.
    Convert a table to df
    Convert gi2operon to df
    -----------------------
    Usage:
    import rotifer.neighborhood.neighborhood as nh
    ipt = nh.ipt2table(str)
    df = ipt.df
    -----------------------
    Other information:
    It is possible to read a table or gi2operons file
    example:
    ipt = nh.ipt2table(open('gi2operon.file').read().splitlines())
    ipt = nh.ipt2table(open('table.file').read().splitlines())
    '''

    def __init__(self, ipt, verbose = False):

        self.ipt = ipt
        self.verbose = verbose
        self.df = self._prepareDF()
    def _prepareDF(self):

        if 'ORGANISM' in self.ipt[0]:
            self.input_format_type = 'gi2operon'
        else:
            self.input_format_type = 'table'

        if self.verbose:
            cf.vmsg('Input format: {0}'.format(self.input_format_type))


        if self.input_format_type == 'table':
            df = self._table2df()
            return df

        if self.input_format_type == 'gi2operon':
            df = self._gi2table()
            return df

    def _table2df(self):
        s = time.time()

        if self.verbose:
            cf.vmsg('Loading table')

        df = pd.read_csv(StringIO("\n".join(self.ipt)), sep = "\t")
        # SLOW!
        # to_int = ['plen', 'start', 'end', 'block_id', 'strand']
        # for col in to_int:
        #     if col == 'start' or col == 'end':
        #         df.astype({col:str})
        #         df['col2'] = np.where(np.core.defchararray.startswith(df[col].astype(str), '>'), 'HEHE', 'HAHA' )
        #         print(df[[col,'col2']].astype(str).values)
        #         # sys.stderr.write('\n'.join('\t'.join(df[[col,'col2']].astype(str).values)))#.values.lstrip('>'), sys.stderr.write(df[col]) )
        # sys.exit()
        # for col in to_int:
        #     if self.verbose:
        #         sys.stderr.write('## converting {0}\n'.format(col))
        #     # df[['a', 'b']] = df[['a','b']].fillna(value=0)
        #     df = df.astype({col: str})
        #     df[[col]] = df[[col]].fillna(0)
        #     try:
        #         df[col] = df[col].map(lambda x: x.lstrip('>'))
        #     except: sys.stderr.write(col)
        #
        #     df[col] = df.astype({col:int})

        df = df.fillna('.')

        try:
            if '' in df['plen'].values:
                df.loc[df[df['plen'] == ''].index, 'plen'] = '.'
        except: pass

        try:
            if '' in df['pid'].values:
                df.loc[df[df['pid'] == ''].index, 'pid'] = '.'
        except: pass

        df = df.replace('', '.')
        df = df.fillna('.')

        if self.verbose:
            cf.vmsg('Loaded DataFrame')

        e = time.time()
        if self.verbose:
            cf.vmsg('Loaded Time: {0:.2f}s'.format(e-s))

        return df


    def _gi2table(self):
        s = time.time()
        block_id = 0
        header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene product organism classification'.split()
        check = 0
        df = pd.DataFrame()
        for line in self.ipt:
            if 'ORGANISM' in line:
                if df.empty and check == 1:
                    df = pd.DataFrame(np.column_stack([nucleotide_ls,
                                       start_ls,
                                       end_ls,
                                       strand_ls,
                                       block_id_ls,
                                       query_ls,
                                       pid_ls,
                                       cds_type_ls,
                                       length_ls,
                                       locus_ls,
                                       chr_type_ls,
                                       asm_ls,
                                       gene_ls,
                                       product_ls,
                                       organism_ls,
                                       classification_ls
                                       ]), columns = header)

                if not df.empty and check == 1:
                    _ = pd.DataFrame(np.column_stack([nucleotide_ls,
                                       start_ls,
                                       end_ls,
                                       strand_ls,
                                       block_id_ls,
                                       query_ls,
                                       pid_ls,
                                       cds_type_ls,
                                       length_ls,
                                       locus_ls,
                                       chr_type_ls,
                                       asm_ls,
                                       gene_ls,
                                       product_ls,
                                       organism_ls,
                                       classification_ls
                                       ]), columns = header)
                    df = pd.concat([df,_])

                block_id +=1
                # Init lists
                chr_type_ls = []
                asm_ls = []
                classification_ls = []
                nucleotide_ls = []
                query_ls = []
                start_ls = []
                end_ls = []
                strand_ls = []
                length_ls = []
                pid_ls = []
                cds_type_ls = []
                gene_ls = []
                locus_ls = []
                gi_ls = []
                product_ls = []
                block_id_ls = []
                organism_ls = []

                chr_type = '.'
                asm = '.'
                classification = '.'
                organism = line.split('ORGANISM')[1].split(' accession no')[0].strip()
                nucleotide = line.split('ORGANISM')[1].split(' accession no is')[1].split('Protein')[0].strip()

                if '|' in line:
                    other_information = line.split('|')
                    for ele in other_information:
                        if 'seq_type' in ele:
                            chr_type = ele.split(':')[1].strip()
                        if 'assembly' in ele:
                            asm = ele.split(':')[1].strip()
                        if 'classification' in ele:
                            classification = ele.split(':')[1].strip()

            if 'ORGANISM' not in line and 'cds' not in line and '---------------------' not in line:
                line = [x for x in line.split(' ') if x != '']
                query, start_end, strand, length, pid,cds_type, gene, locus, gi, *product = line
                start, end = start_end.split('..')
                query = 1 if query == '-->' else 0
                strand = 1 if strand == '+' else -1

                nucleotide_ls.append(nucleotide)
                start_ls.append(start)
                end_ls.append(end)
                strand_ls.append(strand)
                block_id_ls.append(block_id)
                query_ls.append(query)
                pid_ls.append(pid)
                cds_type_ls.append(cds_type)
                length_ls.append(length)
                locus_ls.append(locus)
                chr_type_ls.append(chr_type)
                asm_ls.append(asm)
                gene_ls.append(gene)
                product_ls.append(' '.join(product))
                organism_ls.append(organism)
                classification_ls.append(classification)
                check = 1

        if not df.empty and check == 1:
            _ = pd.DataFrame(np.column_stack([nucleotide_ls,
                               start_ls,
                               end_ls,
                               strand_ls,
                               block_id_ls,
                               query_ls,
                               pid_ls,
                               cds_type_ls,
                               length_ls,
                               locus_ls,
                               chr_type_ls,
                               asm_ls,
                               gene_ls,
                               product_ls,
                               organism_ls,
                               classification_ls
                               ]), columns = header)
            df = pd.concat([df,_])

        if df.empty and check == 1:
            _ = pd.DataFrame(np.column_stack([nucleotide_ls,
                               start_ls,
                               end_ls,
                               strand_ls,
                               block_id_ls,
                               query_ls,
                               pid_ls,
                               cds_type_ls,
                               length_ls,
                               locus_ls,
                               chr_type_ls,
                               asm_ls,
                               gene_ls,
                               product_ls,
                               organism_ls,
                               classification_ls
                               ]), columns = header)
            df = _
        if self.verbose:
            cf.vmsg('Loaded DataFrame')
        e = time.time()
        if self.verbose:
            cf.vmsg('Loaded Time: {0}'.format(str(e-s)))

        return df.drop_duplicates().reset_index(drop = True)

class writer:
    '''
    Print a pandas dataframe in gi2operon format or a table format
    Parameters:
    df = Pandas DataFrame
    of = Output format (table or gi2operon) [Default: table]
    other_info = Other information in gi2operon format
    '''
    def __init__(self,df, of = 'table', other_info='', new_info_ls = '', pad = False,
                 exclude_col = '', include_col = '', verbose = False,
                 print_header = False, **kwargs):
        self.verbose = verbose
        # sys.stderr.write('\n'.join([str(x) for x in df.dtypes]))

        self.print_header = print_header
        try:
            self.kwargs = kwargs['kwargs']
        except:
            self.kwargs = ''

        # sys.stderr.write('haha')
        if verbose:
            cf.vmsg('Preparing to print')
        if of == 'gi2operon' or of == 'gi2operons':
            self.gi2operon(df, other_info, new_info_ls, pad,
                           exclude_col, include_col)

        if of == 'table':
            self.table(df,exclude_col, include_col)

        if of == 'compact':
            self.compact(df,other_info, new_info_ls, exclude_col, include_col,
                         pad)

    def compact(self, df, other_info, new_info_ls, exclude_col,
                include_col,
                pad = False):

        if 'arch' in self.kwargs.keys():
            col_arch = self.kwargs['arch']

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

        new_info_ls = [x for x in new_info_ls if x not in exclude_col]
        if new_info_ls:
            for k in new_info_ls:
                dc[k] = self._get_block_id

        g = df.groupby('block_id').agg(dc).reset_index()

        if pad:
            g['block_id'] = g['block_id'].astype(str)
            g['left'] = g['arch_compact'].str.split('*').str[0]

            g['right'] = g['arch_compact'].map(lambda x: '*'+''.join(x.split('*')[1:] if len(x.split('*')) >1 else '*'))
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

            g['right'] = g['arch_compact'].map(lambda x: '*'+''.join(x.split('*')[1:] if len(x.split('*')) >1 else '*'))
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


    def gi2operon(self,df, other_info = '', new_info_ls= '', pad = False,
                  exclude_col = '',
                  include_col = ''):

        df = df.reset_index(drop = True)

        df = df.astype(str)

        df["block_id"] = pd.to_numeric(df["block_id"])

        to_str = ['plen', 'start', 'end', 'strand', 'query']
        for col in to_str:
            try:
                df = df.astype({col: str})
            except: pass

        if exclude_col:
            if 'gi' in exclude_col:
                exclude_col.append('modified')
            exclude_col = [x for x in exclude_col if x in df.columns and x not in include_col]
            df = df.drop(columns = exclude_col)

        if include_col:
            include_col = [x for x in include_col if x in df.columns]
            df = df[include_col]

        if self.verbose:
            cf.vmsg('Preparing groups')

        # if pad:
        df['plen'].astype(str)
        if '' in df['plen'].values:
            df.loc[df[df['plen'] == ''].index, 'plen'] = '.'

        if '' in df['pid'].values:
            df.loc[df[df['pid'] == ''].index, 'pid'] = '.'
        col_len = self.collen(df, new_info_ls)

        # if no pad
        # TODO

        header = ['.']
        if 'start' not in exclude_col or 'end' not in exclude_col:
            header.append('cds')

        if 'strand' not in exclude_col:
            header.append('dir')

        if 'len' not in exclude_col:
            header.append('len')

        if 'pid' not in exclude_col:
            header.append('pid')

        if 'type' not in exclude_col:
            header.append('type')

        if 'gene' not in exclude_col:
            header.append('gene')

        if 'locus' not in exclude_col:
            header.append('locus')

        if 'gi' not in exclude_col:
            header.append('gi')

        if 'product' not in exclude_col:
            header.append('product')

        if new_info_ls:
            header.extend(new_info_ls)

        header = [header[int(i)].ljust(int(col_len[int(i)])) for i in range(0,len(header), 1)]

        # Preparing
        df['dot'] = np.where(df['query'] == '1', '-->', '.')
        df['dir'] = np.where(df['strand'] == '1', '+', '-')
        df['cds'] = df['start'].astype(str) + '..' + df['end'].astype(str)
        df['len'] =df['plen']
        try:
            df['gi'] = df['modified']
        except:
            df['gi'] = ['.'] * df.shape[0]

        cols_toprint = ['dot']

        if 'start' not in exclude_col or 'end' not in exclude_col:
            cols_toprint.append('cds')
        if 'strand' not in exclude_col:
            cols_toprint.append('dir')
        if 'len' not in exclude_col:
            cols_toprint.append('len')
        if 'pid' not in exclude_col:
            cols_toprint.append('pid')
        if 'type' not in exclude_col:
            cols_toprint.append('type')
        if 'gene' not in exclude_col:
            cols_toprint.append('gene')
        if 'locus' not in exclude_col:
            cols_toprint.append('locus')
        if 'gi' not in exclude_col:
            cols_toprint.append('gi')
        if 'product' not in exclude_col:
            cols_toprint.append('product')

        if new_info_ls:
            for new_col in new_info_ls:
                cols_toprint.append(new_col)

        # Merge information
        if self.verbose:
            cf.vmsg('Combining information')
        df['print_here'] = df.apply(lambda x: self.combine_cols_for_gi2operon(x, cols_toprint,col_len ), 1)

        if self.verbose:
            cf.vmsg('All information collected')
        # Group by Block_id

        #agg

        dc_group = {'print_here':'sum',
                'pid': self.get_protein_in_block,
                'query':self.get_protein_in_block,
                'nucleotide': self.get_protein_in_block,
                'organism': self.get_protein_in_block}

        if other_info:
            for e in other_info:
                if e != 'block_id':
                    dc_group[e] = self._get_block_id

        g = df.groupby('block_id').agg(dc_group).reset_index()

        total_groups = g.shape[0]
        x = 0

        for i, sub_df in g.iterrows():
            if self.verbose:
                cf.vmsg('Remaining groups: {0}'.format(total_groups-x))
            x +=1

            dc_sub_df = {'query': sub_df['query'],
                              'pid': sub_df['pid'],
                              'nucleotide':sub_df['nucleotide'],
                              'organism': sub_df['organism']}
            if other_info:
                for e in other_info:
                    dc_sub_df[e] = sub_df[e]
            _ = pd.DataFrame(dc_sub_df)

            try:
                protein = _[_['query'] == '1']['pid'].values[0]
            except: protein = '.'
            nucleotide = _.nucleotide.values[0]

            max_header = 'ORGANISM ' + _.organism.values[0] + ' accession no is ' + nucleotide + ' Protein is '+ protein
            if other_info:
                for e in other_info:
                    max_header += ' | '+ e + ':'+''.join(str(_[e].values[0]))

            print(max_header)
            print('  '.join(header))
            print(sub_df['print_here'].rstrip('\n'))
            print("---------------------------------------")

    def get_protein_in_block(self, row):
        # print(row.values)
        return list(row.values)


    def combine_cols_for_gi2operon(self, row, cols, col_len):
        ls = []
        for x in range(len(cols)):
            ls.append(row[cols[x]])

        return '  '.join([ls[i].ljust(int(col_len[i])) for i in range(len(ls))]) +'\n'

    def table(self, df, exclude_col = '', include_col = ''):

        if exclude_col:
            if 'gi' in exclude_col:
                exclude_col.append('modified')
            exclude_col = [x for x in exclude_col if x in df.columns and x not in include_col]
            df = df.drop(columns = exclude_col)

        if include_col:
            include_col = [x for x in include_col if x in df.columns]
            df = df[include_col]

        if self.print_header:
            header = df.columns
            print('\t'.join(header))

        else:
            pass
        for i,row in df.iterrows():
            print('\t'.join(row.map(str).values))

    #######################
    #### Column length ####
    #######################

    def collen(self,sub_df, new_info_ls =''):
        '''
        Pad for gi2operon
        '''
        collen_len = []
        try: sub_df['start'].astype(str)
        except: pass
        try: sub_df['end'].astype(str)
        except: pass
        try: sub_df['plen'].astype(str)
        except: pass

        q = '-->'
        collen_len.append(len(q)) # query mark

        try:
            sub_df['cds_loc'] = sub_df['start'].astype(str)+'..'+sub_df['end'].astype(str)
        except: pass

        try:
            len_cds = sub_df['cds_loc'].str.len().max()
            collen_len.append(len_cds)
        except:
            pass

        direction = len('dir')
        collen_len.append(direction)

        try:
            len_plen = sub_df['plen'].astype(str).str.len().max() if sub_df['plen'].astype(str).str.len().max() >= len('len') else len('len')
            collen_len.append(len_plen)
        except: pass

        try:
            len_pid = sub_df['pid'].str.len().max() if sub_df['pid'].str.len().max() >= len('pid') else len('pid')
            collen_len.append(len_pid)
        except: pass

        try:
            len_type = sub_df['type'].str.len().max() if sub_df['type'].str.len().max() >= len('type') else len('type')
            collen_len.append(len_type)
        except: pass

        try:
            len_gene = sub_df['gene'].str.len().max() if sub_df['gene'].str.len().max() >= len('gene') else len('gene')
            collen_len.append(len_gene)
        except:
            pass

        try:
            len_locus = sub_df['locus'].str.len().max() if sub_df['locus'].str.len().max() >= len('locus') else len('locus')
            collen_len.append(len_locus)
        except: pass

        try:
            try:
                len_modified = sub_df['modified'].str.len().max() if sub_df['modified'].str.len().max() >= len('gi') else len('gi')
                collen_len.append(len_modified)

            except:
                len_modified = 1
                collen_len.append(len_modified)
        except: pass

        try:
            len_product = sub_df['product'].str.len().max() if sub_df['product'].str.len().max() >= len('product') else len('product')
            collen_len.append(len_product)
        except: pass

        if new_info_ls:
            for col in new_info_ls:
                new_col = sub_df[col].str.len().max() if sub_df[col].str.len().max() >= len(col) else len(col)
                collen_len.append(new_col)

        return collen_len

    ############################
    #### <\ Column length > ####
    ############################

    ###########################
    #### < _get_block_id > ####
    ###########################

    def _get_block_id(self, series):
        query = series.unique()
        if len(query) ==1:
            return query
        else:
            query = '; '.join(list(series.unique()))
            return query

    ##################################

    def old_gi2operon(self):
        df = df.reset_index(drop = True)

        df = df.astype(str)

        df["block_id"] = pd.to_numeric(df["block_id"])

        to_str = ['plen', 'start', 'end', 'strand', 'query']
        for col in to_str:
            try:
                df = df.astype({col: str})
            except: pass

        if exclude_col:
            if 'gi' in exclude_col:
                exclude_col.append('modified')
            exclude_col = [x for x in exclude_col if x in df.columns and x not in include_col]
            df = df.drop(columns = exclude_col)

        if include_col:
            include_col = [x for x in include_col if x in df.columns]
            df = df[include_col]

        if self.verbose:
            sys.stderr.write('Preparing groups')
        g = df.groupby('block_id')

        if pad:
            df['plen'].astype(str)
            if '' in df['plen'].values:
                df.loc[df[df['plen'] == ''].index, 'plen'] = '.'

            if '' in df['pid'].values:
                df.loc[df[df['pid'] == ''].index, 'pid'] = '.'
            col_len = self.collen(df, new_info_ls)

        for i, sub_df in g:
            sub_df['block_id'] = sub_df['block_id'].astype(str)

            if not pad:
                sub_df['plen'].astype(str)
                if '' in sub_df['plen'].values:
                    sub_df.loc[sub_df[sub_df['plen'] == ''].index, 'plen'] = '.'

                if '' in sub_df['pid'].values:
                    sub_df.loc[sub_df[sub_df['pid'] == ''].index, 'pid'] = '.'
                col_len = self.collen(sub_df, new_info_ls)

            try:
                protein = sub_df[sub_df['query'] == '1']['pid'].values[0]
            except: protein = '.'
            nucleotide = sub_df.nucleotide.values[0]

            max_header = 'ORGANISM ' + sub_df.organism.values[0] + ' accession no is ' + nucleotide + ' Protein is '+ protein
            if other_info:
                for e in other_info:
                    max_header += ' | '+ e + ':'+''.join(str(sub_df[e].values[0]))

            print(max_header)

            #sub header
            # exlude from header
            header = ['.']
            if 'start' not in exclude_col or 'end' not in exclude_col:
                header.append('cds')

            if 'strand' not in exclude_col:
                header.append('dir')

            if 'len' not in exclude_col:
                header.append('len')

            if 'pid' not in exclude_col:
                header.append('pid')

            if 'type' not in exclude_col:
                header.append('type')

            if 'gene' not in exclude_col:
                header.append('gene')

            if 'locus' not in exclude_col:
                header.append('locus')

            if 'gi' not in exclude_col:
                header.append('gi')

            if 'product' not in exclude_col:
                header.append('product')

            if new_info_ls:
                header.extend(new_info_ls)

            header = [header[int(i)].ljust(int(col_len[int(i)])) for i in range(0,len(header), 1)]
            print('  '.join(header))

            # sub_df['dot'] = np.where(sub_df['query'] == '1', '-->', '.')
            sub_df['.'] = np.where(sub_df['query'] == '1', '-->', '.')
            sub_df['dir'] = np.where(sub_df['strand'] == '1', '+', '-')
            sub_df['cds'] = sub_df['start'].astype(str) + '..' + sub_df['end'].astype(str)
            # sub_df['protlen'] =sub_df['plen']
            sub_df['len'] =sub_df['plen']
            sub_df['gi'] = sub_df.shape[0]*['.']

            for i, row in sub_df.iterrows():
                toprint = [row['.']]
                if 'start' not in exclude_col or 'end' not in exclude_col:
                    toprint.append(str(row.cds))
                if 'strand' not in exclude_col:
                    toprint.append(row['dir'])
                if 'len' not in exclude_col:
                    toprint.append(row['len'])
                if 'pid' not in exclude_col:
                    toprint.append(row['pid'])
                if 'type' not in exclude_col:
                    toprint.append(row['type'])
                if 'gene' not in exclude_col:
                    toprint.append(row['gene'])
                if 'locus' not in exclude_col:
                    toprint.append(row['locus'])
                if 'gi' not in exclude_col:
                    toprint.append(row['gi'])
                if 'product' not in exclude_col:
                    toprint.append(row['product'])

                if new_info_ls:
                    for new_col in new_info_ls:
                        toprint.append(row[new_col])
                toprint = [toprint[i].ljust(int(col_len[i])) for i in range(len(toprint))]
                print('  '.join(toprint))
            print("---------------------------------------")
