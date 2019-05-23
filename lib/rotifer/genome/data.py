#!/usr/bin/env python3

# data is a child of Pandas DataFrame manipulate DataFrame here
#

import pandas as pd
import rotifer.core.log as rlog

class NeighborhoodDF(pd.DataFrame):
    # mandatory columns

    _metadata = ['filterby',
                 'verbose',
                 'log_file']

    def __init__(self,
                 *args, **kwargs):
        self.filterby = kwargs.pop('filterby', None)
        self.verbose = kwargs.pop('version', None)
        self.log_file = kwargs.pop('log_file', None)
        super(NeighborhoodDF, self).__init__(*args, **kwargs)
        if self.filterby:
            # do stuf
            pass

    # Pandas contructor

    @property
    def _constructor(self):
        return NeighborhoodDF

    def writer(self,
               **kwargs):

        '''
        kwargs:
            - of/output_format: Output format (default: table)
                                - table
                                - pretty
                                - gi2operon
                                - gi2op
                                - gi2operons
        '''

        self._of = 'table'
        if 'of' in kwargs.keys() or 'output_format' in kwargs.keys():
            try:
                self._of = kwargs['of'].lower()
            except:
                self._of = kwargs['output_format'].lower()

        if self._of.lower() in ['gi2operon', 'gi2operons', 'gi2op']:
            self._gi2operon(**kwargs)

        if self._of == 'table':
            self._table()

        if self._of == 'compact':
            self._compact()

    def _gi2operon(self, **kwargs):

        self._exclude_col = []
        self._include_col = []
        self._new_info_ls = []

        # Put df
        df = self

        df = df.reset_index(drop = True)

        df = df.astype(str)

        # Convert block id to number
        df["block_id"] = pd.to_numeric(df["block_id"])

        to_str = ['plen', 'start', 'end', 'strand', 'query']
        for col in to_str:
            try:
                df = df.astype({col: str})
            except: pass

        if self._exclude_col:
            if 'gi' in self._exclude_col:
                exclude_col.append('modified')

            self._exclude_col = [x for x in self._exclude_col if x in df.columns and x not in self._include_col]

            df = df.drop(columns = self._exclude_col)

        if self._include_col:
            include = []
            for x in self._include_col:
                if x in df.columns:
                    include.append(x)

                else:
                    rlog.log({
                         2: f"Column named [{columns}] was not included the DataFrame"},
                         level = self._verbose,
                         name = __name__,
                         log_file = self._log_file)

            df = df[include]

        try:
            df['plen'] = df['plen'].map(lambda x: '.' if x == '' else x)

        except:
            pass

        try:
            df['pid'] = df['pid'].map(lambda x: '.' if x == '' else x)

        except:
            pass

        col_len = self._collen(df)

    def _collen(self, df):
        collen_len = []
        try: self.df['start'].astype(str)
        except: pass
        try: df['end'].astype(str)
        except: pass
        try: df['plen'].astype(str)
        except: pass

#        self._new_info_ls

        q = '-->'

        collen_len.append(len(q)) # query mark

        try:
            df['cds_loc'] = sub_df['start'].astype(str)+'..'+sub_df['end'].astype(str)
        except: pass

        try:
            len_cds = df['cds_loc'].str.len().max()
            collen_len.append(len_cds)

        except:
            pass

        direction = len('dir')
        collen_len.append(direction)

        try:
            len_plen = df['plen'].astype(str).str.len().max() if df['plen'].astype(str).str.len().max() >= len('len') else len('len')
            collen_len.append(len_plen)

        except: pass

        try:
            len_pid = df['pid'].str.len().max() if df['pid'].str.len().max() >= len('pid') else len('pid')
            collen_len.append(len_pid)
        except: pass

        try:
            len_type = df['type'].str.len().max() if df['type'].str.len().max() >= len('type') else len('type')
            collen_len.append(len_type)

        except: pass

        try:
            len_gene = df['gene'].str.len().max() if df['gene'].str.len().max() >= len('gene') else len('gene')
            collen_len.append(len_gene)

        except:
            pass

        try:
            len_locus = df['locus'].str.len().max() if df['locus'].str.len().max() >= len('locus') else len('locus')
            collen_len.append(len_locus)
        except: pass

        try:
            try:
                len_modified = df['modified'].str.len().max() if df['modified'].str.len().max() >= len('gi') else len('gi')
                collen_len.append(len_modified)

            except:
                len_modified = 1
                collen_len.append(len_modified)
        except: pass

        try:
            len_product = df['product'].str.len().max() if df['product'].str.len().max() >= len('product') else len('product')
            collen_len.append(len_product)
        except: pass

        if self._new_info_ls:
            for col in self._new_info_ls:
                new_col = df[col].str.len().max() if df[col].str.len().max() >= len(col) else len(col)
                collen_len.append(new_col)

        return collen_len

    def _table(self, header = False):

        pass

    def _compact(self):
        pass

    def concat(self, other):
        return pd.concat([self, other])






    # FILTERS
