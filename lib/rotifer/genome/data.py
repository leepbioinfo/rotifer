#!/usr/bin/env python3

# data is a child of Pandas DataFrame manipulate DataFrame here
#
import itertools
import numpy as np
import pandas as pd
import rotifer.core.log as rlog


def pandas_print_everything():
    """
    A small function to allow pandas to print not only a few columns and rows.
    Take care when use it, big dataframes could take a long time to print.
    """
    pd.options.display.max_rows = 10000000000
    pd.options.display.max_columns = 1000000000
    pd.options.display.max_colwidth = 1000000000
    pd.options.display.width = 100000
    return print('Pandas printing the full dataframe')


def not_kwargs(dict_args, key, value):
    """ Small function to check if a element is not in a dictionary.
    It is useful to set a default value to functions that exceeds the
    number of 5 arguments
    """
    if key in dict_args.keys():
        return dict_args[key]
    return value

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

    def _table(self):
        pass

    def _compact(self):
        pass



###### New functions added 12/08/2019

    def group_strand(self):
        """
        A function to add the column strand_groups in the acc2operon dataframe.
        """
        dataframe = self.copy()
        dataframe['strand_group'] = np.cumsum(
            np.where(dataframe.block_id == dataframe.block_id.shift(1),
                     np.where(
                         dataframe.strand.shift(1) == dataframe.strand, 0, 1
                     ), 1)
        )
        return dataframe

    def same_strand(self,
                    query='query'):
        """
        Using the strand group and a query it will filter the self with
        strands that present one  query at least.
        """
        dataframe = self.group_strand().copy()
        return dataframe[dataframe.strand_group.isin(
            dataframe.groupby('strand_group')[query].sum().where(
                lambda x: x > 0
            ).dropna().index)]

    def value_counts_full(self,
                          column):
        """
        Calculate absolute and relative frequencies for a given column.
        If dropna is True, it will ignore the empty values to count
        the relative frequency.
        """
        x = self[self['type']=="CDS"]
        lpid = x['pid'].nunique()
        lbid = x['block_id'].nunique()
        g_df = x.groupby(column).agg({'block_id':'nunique','pid':'count'}).sort_values('block_id', ascending=False).reset_index()
        g_df['relative_block_id'] = g_df.block_id/lbid *100
        g_df['relative_pid'] = g_df.pid/lpid * 100
        g_df.rename({'block_id':'block_id_count', 'pid':'pid_count'}, axis=1)
        print(f'Total of block_id in selection:{lbid}')
        print(f'Total of pid in selection :{lpid}')

        return g_df

    def series_to_compact_frequency(self, concat_string=": ", normalize=False):
        """
         *Series method*
        Aggregate the value count of a series in a single string with
        the count inside a parentheses
        """
        if normalize:
            count_dict = dict(
                self.value_counts(normalize=normalize).map('{:,.2f}'.format)
            )
        else:
            count_dict = dict(self.value_counts(normalize=normalize))
        count_list = []
        for each_series in count_dict:
            count_list.append(f'{each_series}({count_dict[each_series]})')
        return concat_string.join(count_list)

    def select_neighbors(self,
                         column='query',
                         value=1,
                         new_query=True,
                         **kwargs):
        """
        Slice the self in x neighbors of a given value query.
        self : The Dataframe to be sliced.
        Column : The column that will be filtered.
        value: The value to be filtered.
        after and before: The numbers of neighbors filtered.
        strand:Filter for neighbors that are at the same strand of the query
        stats: If set it will give the stats of the column
        Examples:
        get_x_neighbors(DF_ATPases, column='Pfam', value='AAA_9',before=3,after=3).
        In this example the Dataframe DF_ATPases will be filtered by only blocks
        that present the AAA_9 ATPase and will return 3 before and 3 after neighbors
        get_x_neighbors(DF_ATPases, column='cluster', value='sc_29',
        before=3,after=3, stats='pfam).
        The output of this example will be the stats of pfam column
        in the neighborhood of cluster sc_29.
        """

        after = not_kwargs(kwargs, 'after', 1)
        before = not_kwargs(kwargs, 'before', 1)
        strand = not_kwargs(kwargs, 'strand', False)
        stats = not_kwargs(kwargs, 'stats', False)



        query_df = self[self[column] == value].copy()
        full_query_df = \
                self[self.block_id.isin(query_df.block_id)].copy()
        min_max = full_query_df.reset_index().groupby(
            'block_id'
        ).agg({'index': ['min', 'max']})
        min_max.columns = min_max.columns.droplevel()
        query_min_max = \
                query_df.reset_index()[['block_id', 'index']].set_index('block_id')
        query_min_max = query_min_max.join(min_max)
        query_min_max['before_limit'] = 1 + (
            np.where(query_min_max['index'] + before >= query_min_max['max'],
                     query_min_max['max'],
                     query_min_max['index'] + before)
        )
        query_min_max['after_limit'] = np.where(
            query_min_max['index'] -after <= query_min_max['min'],
            query_min_max['min'],
            query_min_max['index'] - after
        )
        query_min_max['Diference'] = \
                query_min_max.before_limit - query_min_max.after_limit
        index_list = list(
            itertools.chain.from_iterable(
                [list(range(x[0], x[1])) for x in list(
                    zip(query_min_max.after_limit, query_min_max.before_limit))]))
        selected_df = self.iloc[index_list].copy()
        selected_df.loc[:, 'new_query'] = 0
        selected_df.loc[selected_df[column] == value, 'new_query'] = 1
        if strand:
            if stats:
                return selected_df.same_strand(query='new_query').value_counts_full(stats)
            if new_query:
                return selected_df.same_strand(query='new_query')
            return selected_df.same_strand(query='new_query').drop(
                'new_query', axis=1)

        if stats:
            return selected_df.value_counts_full(stats)
        if new_query:
            return selected_df
        return selected_df.drop('new_query', axis=1)

    def compact_neighborhood(self,
                             annotation='pfam',
                             query='query',
                             strand=False,
                             fill_unk=False):
        """
        A tabular and compact representation of  selected neighborhoods block.
        Given a acc2operons Dataframe it will print the block ID as first column
        and a compact form of all locus within the block_id as compact form using
        a given
        column to annotate it.
        If stand is True, it will display only the loci that are at the same
        strand of the query, and also the ID displayed is a new calculated Strand ID.
        """

        def compacted(dataframe, annotation=annotation):
            """
            Take a groupby slice of the dataframe and make a new collum
            containing  a string of  annotation  and strand direction.
            """

            dataframe[annotation].fillna('?', inplace=True)
            dataframe[annotation] = np.where(
                dataframe['type'] != 'CDS',
                dataframe['type'],
                dataframe[annotation]
            )
            if fill_unk:
                dataframe[annotation] = np.where(
                    dataframe[annotation] == '?',
                    dataframe[fill_unk],
                    dataframe[annotation])

            dataframe[annotation] = dataframe[annotation].astype(str)
            return np.where(
                dataframe['strand_group'] == dataframe['strand_group'].shift(-1),
                np.where(
                    dataframe[query] == 0,
                    np.where(
                        dataframe.strand == 1,
                        dataframe[annotation]+'->',
                        '<-'+dataframe[annotation]),
                    np.where(
                        dataframe.strand == 1,
                        '*' + dataframe[annotation]+'->',
                        '<-'+dataframe[annotation]+'*')
                ),
                np.where(
                    dataframe[query] == 0,
                    np.where(
                        dataframe.strand == 1,
                        dataframe[annotation]+'->|',
                        '<-'+dataframe[annotation] + '|'),
                    np.where(
                        dataframe.strand == 1,
                        '*' + dataframe[annotation]+'->|',
                        '<-'+dataframe[annotation]+ '*|'),
                )
            )

        def main_cnei(df_to_compact):
            """
            Main function to make the compact_neighborhood.
            It was designed to be used in a groupby.
            """
            #Checking the orientation of the first query, and if it is
            # in the reverse strand, change the orientation. This way all
            # the blocks will be displayed at the same orientation.
            #Somehow if you are trying to groupby one block without a query
            # an error will be raised.
            try:
                if df_to_compact[df_to_compact[query] == 1].loc[:, 'strand'].iloc[0] == 1:
                    _ = df_to_compact.copy()
                else:
                    _ = df_to_compact.copy().iloc[::-1]
                    _.strand = _.strand *- 1

            except IndexError:
                raise NameError('Query_Missing')

            #Fill the no annotated with an ? marker or a column in the Dataframe
            # that the users can select by the option fill_unk=column name.

            _['compact'] = compacted(_)
            _['pid_compact'] = compacted(_, annotation='pid')


            return _.groupby(
                group
            ).agg({'compact': lambda x: "".join(list(x))[:-1],
                   'block_id':'first',
                   'lineage':'first',
                   'organism':'first',
                   'seq_type':'first',
                   'classification':'first',
                   'pid_compact': lambda x: "".join(list(x))[:-1]
                  })


        # Main function, select if it will be grouped by block_id or strand,
        # and then make the groupby.
        group = 'block_id'
        if strand:
            group = 'strand_group'
            _1 = self.same_strand(query=query).copy()
        else:
            _1 = self.group_strand().copy()


        return _1.groupby(group).apply(main_cnei).reset_index(
            drop = True).set_index('block_id')


    def full_neighborhood(self, column='cluster', value='value', stats=False, new_query=True):
        """
        Function that filter a given value in a collumn and fetch
        all the loci within the block_id
        """
        selected_df = self[self.block_id.isin(
            self[self[column] == value].block_id
        )].copy()
        if stats:
            return  selected_df.value_counts_full(stats)

        selected_df.loc[:, 'new_query'] = 0
        selected_df.loc[selected_df[column] == value, 'new_query'] = 1
        selected_df = selected_df[selected_df.block_id.isin(
            selected_df[selected_df[column] == value].block_id)]
        if new_query:
            return selected_df
        return selected_df.drop('new_query', axis=1)




        # FILTERS
