#!/usr/bin/env python3

# Data here is a child of Pandas DataFrame

import itertools
import numpy as np
import pandas as pd
import rotifer.core.log as rlog
from rotifer.core.functions import not_kwargs
from rotifer.taxonomy.utils import lineage as rtlineage

class NeighborhoodDF(pd.DataFrame):
    # mandatory columns????
    _metadata = ['filterby',
                 'verbose',
                 'log_file']

    def __init__(self, *args, **kwargs):
        self.filterby  = kwargs.pop('filterby', None)
        self.verbose   = kwargs.pop('verbose', None)
        self.log_file  = kwargs.pop('log_file', None)
        update_lineage = kwargs.pop('update_lineage',False)
        preferred_taxa = kwargs.pop('preferred_taxa',None)
        super(NeighborhoodDF, self).__init__(*args, **kwargs)

        # Update or add lineage column
        cols = self.columns.to_list()
        if isinstance(update_lineage,str) & (update_lineage in cols):
            if "lineage" in cols:
                self.lineage = rtlineage(self[update_lineage], preferred_taxa=preferred_taxa)
            else:
                taxonomy = [ i for i in range(0,len(cols)) if cols[i] == update_lineage ][0]
                self.insert(taxonomy, "lineage", rtlineage(self[update_lineage], preferred_taxa=preferred_taxa))

    # Pandas contructor
    @property
    def _constructor(self):
        return NeighborhoodDF

    def writer(self, output_format='table', **kwargs):
        '''
        kwargs:
          of | output_format: Output format (default: table)
           Supported output formats:
            - table   : pandas supported tabular output
            - compact : single row description of gene neighborhoods
            - gi2operons | gi2operon | gi2op : old TASS text view format
        '''

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
        df = self.copy()
        df = df.reset_index(drop = True)
        df = df.astype(str)

        # Convert block id to numeric and other columns to string
        df["block_id"] = pd.to_numeric(df["block_id"])
        to_str = ['plen', 'start', 'end', 'strand', 'query']
        for col in to_str:
            try:
                df = df.astype({col: str})
            except: pass
        try:
            df['plen'] = df['plen'].map(lambda x: '.' if x == '' else x)
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
                    rlog.log({2: f"Column named [{columns}] was not included the DataFrame"},
                         level = self._verbose,
                         name = __name__,
                         log_file = self._log_file)
            df = df[include]

        try:
            df['pid'] = df['pid'].map(lambda x: '.' if x == '' else x)

        except:
            pass

        col_len = self._collen(df)

    def _table(self):
        pass

    def _compact(self):
        pass

    def _collen(self, df):
        collen_len = []
        try: self.df['start'].astype(str)
        except: pass
        try: df['end'].astype(str)
        except: pass
        try: df['plen'].astype(str)
        except: pass

        q = '-->'
        collen_len.append(len(q)) # query mark

        try:
            df['cds_loc'] = sub_df['start'].astype(str)+'..'+sub_df['end'].astype(str)
        except: pass

        try:
            len_cds = df['cds_loc'].str.len().max()
            collen_len.append(len_cds)
        except: pass

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

    def value_counts_full(self, column):
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

            #Fill the unannotated with an ? marker or a column in the Dataframe
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
        return _1.groupby(group).apply(main_cnei).reset_index(drop = True).set_index('block_id')


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


    # Find internal_id (iimin,iimax) and feature_order (fomin,fomax) boundaries
    def boundaries(self, nucleotide=None):
        """
        Find maximum and minimum values of block_id, feature_order and internal_id.
        Columns feature_order and internal_id are evaluated per nucleotide entry.

        Returns a Pandas dataframe.
        """
        df = self
        if isinstance(nucleotide,list) or isinstance(nucleotide,pd.Series):
            df = self[self.nucleotide.isin(nucleotide)]
        dflim = df.groupby(['assembly','nucleotide','type']).agg({'feature_order':['min','max']})
        dflim.reset_index(inplace=True)
        dflim.columns = ['assembly','nucleotide','type','fomin','fomax']
        iidlim = df.groupby(['assembly','nucleotide']).agg({'internal_id':['min','max']})
        iidlim.reset_index(inplace=True)
        iidlim.columns = ['assembly','nucleotide','iidmin','iidmax']
        bidlim = df.groupby(['assembly']).agg({'block_id':['min','max']})
        bidlim.reset_index(inplace=True)
        bidlim.columns = ['assembly','bidmin','bidmax']
        dflim = dflim.merge(iidlim, left_on=['assembly','nucleotide'], right_on=['assembly','nucleotide'], how='left')
        dflim = dflim.merge(bidlim, left_on=['assembly'], right_on=['assembly'], how='left')
        return dflim

    def vicinity(self, targets=['query == 1'], before=3, after=3, min_block_distance=0, fttype='same', min_block_id=1):
        """
        Locate genomic regions that contain features selected by some criteria.
        
        Each region will be described by assembly, nucleotide and a pair of minimum and
        maximum internal IDs.

        The user may choose a set of rows (targets) as anchors, whose neighbors will be
        evaluated by user-defined parameters, such as feature type, strand or distance.

        Arguments:
         - targets : (list of) boolean pd.Series or rules to select targets.

            Rules will used identify target rows using Pandas's eval method.

            If all rules return True for a given row, the row is a target and
            its neighbors will be identified.

         - before : keep at most this number of features, of the same type as the target,
                    before each target
         - after  : keep at most this number of features, of the same type as the target,
                    after each target

         - min_block_distance : minimum distance between two consecutive blocks
                                Blocks separated by more features 

         - fttype : how to process feature types of neighbors
                    Supported values:
                     - same : consider only features of the same type as the target
                     - any  : ignore feature type and count all features when
                              setting neighborhood boundaries

         - min_block_id : starting number for block_ids, useful if calling this method
                          eqeuntially through multiple rotifer.genome.data.NeighborhoodDF
                          dataframes

         - circular : whether to go around coordinate 1 in circular genomes (boolen)
        """
        import sys
        import numpy as np
        if min_block_id > 0:
            min_block_id -= 1

        # Build boolean pandas.Series to mark targets
        select = True
        if not isinstance(targets, list):
            targets = [targets]
        for code in targets:
            if isinstance(code,str):
                try:
                    select &= self.eval(code)
                except:
                    print(f'Rule {code} failed', file=sys.stderr)
            elif isinstance(code,pd.Series):
                select &= code

        # Find targets and their neighborhood
        if (not isinstance(select,pd.Series)) or (select.sum() == 0):
            print(f'No anchors to search for neighbors were found! Revise your list of targets!')
            return pd.DataFrame()

        # Make sure data is sorted and has compatible internal ids
        if self.is_fragment.any():
            self.sort_values(['assembly','nucleotide','start','end','internal_id'], inplace=True)
            self.internal_id = list(range(0,len(self)))

        # Initialize dataframe for each region (blocks)
        cols = ['assembly','nucleotide','topology','type']
        dflim = self.boundaries()

        # Process features with type
        if (fttype == 'same'):
            blks = self[select].filter([*cols,'block_id','feature_order','is_fragment'])
            blks.sort_values([*cols,'feature_order'], inplace=True)
            blks['foup']   = blks.feature_order - before
            blks['fodown'] = blks.feature_order + after

            # Identify blocks by merging neighborhoods
            bid = (blks.nucleotide == blks.nucleotide.shift(1))
            bid = (blks.is_fragment & (blks.block_id == blks.block_id.shift(1)))
            bid = bid & (blks.type == blks.type.shift(1))
            bid = bid & ((blks.foup - blks.fodown.shift(1) - 1) <= min_block_distance)
            blks['block_id'] = (~bid).cumsum() + min_block_id
            blks = blks.groupby([*cols,'block_id']).agg({'feature_order':list,'foup':min,'fodown':max,'is_fragment':'all'}).reset_index()
            blks = blks.merge(dflim, left_on=['assembly','nucleotide','type'], right_on=['assembly','nucleotide','type'], how='left')
            blks.rename({'feature_order':'targets'}, axis=1, inplace=True)
            blks['origin'] = 0

            # Initiate analysis of circular replicons
            circular = blks.query('topology == "circular"')
            if len(circular) > 0:
                circular = circular.groupby(['assembly','nucleotide','type'])
                circular = circular.agg({'block_id':['min','max'],'foup':'min','fodown':'max','fomin':'min','fomax':'max'})
                circular.reset_index(inplace=True)
                circular.columns = ['assembly','nucleotide','type','bidmin','bidmax','foup','fodown','fomin','fomax']

                # Merge first and last blocks of circular replicons, if too close
                tooclose = circular.eval(f'(fomax - fodown) + (foup - fomin) <= {min_block_distance}')
                if tooclose.any():
                    bidmaxIndex = blks[blks.block_id.isin(circular[tooclose].bidmax)].index
                    bidminIndex = blks[blks.block_id.isin(circular[tooclose].bidmin)].index
                    #print("First and last are too close:\n"+blks.iloc[bidminIndex.append(bidmaxIndex)].to_string()+"\n")
                    blks.loc[bidminIndex,'foup']     = blks.loc[bidminIndex,'fomin'] # Set foup   = fomin
                    blks.loc[bidmaxIndex,'fodown']   = blks.loc[bidmaxIndex,'fomax'] # Set fodown = fomax
                    blks.loc[bidmaxIndex,'block_id'] = blks.loc[bidminIndex,'block_id'] # Set block_ids to be the same
                    blks.loc[bidminIndex,'origin']   = 1
                    blks.loc[bidmaxIndex,'origin']   = 1
                    circular = circular[~tooclose]

                # Neighborhood runs through the entire replicon
                overrun = (circular.foup < circular.fomin) & (circular.fodown > circular.fomax)
                if overrun.any():
                    overIndex = blks[blks.block_id.isin(circular[overrun].bidmin)].index
                    #print("Overrun:\n"+blks.iloc[overIndex].to_string()+"\n")
                    blks.loc[overIndex,'foup']   = blks.loc[overIndex,'fomin'] # Set foup   = fomin
                    blks.loc[overIndex,'fodown'] = blks.loc[overIndex,'fomax'] # Set fodown = fomax
                    blks.loc[overIndex,'origin'] = 2
                    circular = circular[~overrun]

                # Neighborhood starts before the origin of replication
                startsBefore = (circular.foup < circular.fomin)
                if startsBefore.any():
                    beforeOriginIndex = blks[blks.block_id.isin(circular[startsBefore].bidmin)].index
                    blks.loc[beforeOriginIndex,'origin'] = 1
                    beforeOrigin = blks.iloc[beforeOriginIndex].copy()
                    #print("Starts before origin:\n"+beforeOrigin.to_string()+"\n")
                    beforeOrigin.foup   = beforeOrigin.fomax + beforeOrigin.foup - beforeOrigin.fomin + 1
                    beforeOrigin.fodown = beforeOrigin.fomax
                    blks.loc[beforeOriginIndex,'foup'] = blks.loc[beforeOriginIndex,'fomin'] # Set foup = fomin
                    blks = pd.concat([blks, beforeOrigin]).reset_index(drop=True)
                    circular = circular[~startsBefore]

                # Neighborhood ends after the origin of replication
                endsAfter = (circular.fodown > circular.fomax)
                if endsAfter.any():
                    afterOriginIndex = blks[blks.block_id.isin(circular[endsAfter].bidmax)].index
                    blks.loc[afterOriginIndex,'origin'] = 1
                    afterOrigin = blks.iloc[afterOriginIndex].copy()
                    #print("Ends after origin:\n"+afterOrigin.to_string()+"\n")
                    afterOrigin.foup   = afterOrigin.fomin
                    afterOrigin.fodown = afterOrigin.fomin + afterOrigin.fodown - afterOrigin.fomax - 1
                    blks.loc[afterOriginIndex,'fodown'] = blks.loc[afterOriginIndex,'fomax'] # Set fodown = fomax
                    blks = pd.concat([blks, afterOrigin]).reset_index(drop=True)

            # Truncate feature_id boundaries
            blks.foup   = np.where(blks.foup < blks.fomin, blks.fomin, blks.foup)
            blks.fodown = np.where(blks.fodown > blks.fomax, blks.fomax, blks.fodown)

            # Sort again and add row ID, i.e. within block row number
            blks.sort_values([*cols,'block_id','foup'], ascending=[True,True,True,True,True,False], inplace=True)
            blks.reset_index(inplace=True, drop=True)
            blks['rid'] = blks.index.to_series() - (blks.index.to_series() * (blks.block_id != blks.block_id.shift(1))).cummax()

            # Add first and last internal_id per block
            blks = blks.merge(
                    self[['assembly','nucleotide','type','feature_order','internal_id']],
                    left_on=['assembly','nucleotide','type','foup'],
                    right_on=['assembly','nucleotide','type','feature_order'],
                    how='left').rename({'internal_id':'up'}, axis=1).drop('feature_order', axis=1)
            blks = blks.merge(
                    self[['assembly','nucleotide','type','feature_order','internal_id']],
                    left_on=['assembly','nucleotide','type','fodown'],
                    right_on=['assembly','nucleotide','type','feature_order'],
                    how='left').rename({'internal_id':'down'}, axis=1).drop('feature_order', axis=1)

            # Adjust internal_id boundaries for blocks with all features of the same type as the query
            blks.up   = np.where(blks.foup   == blks.fomin, blks.iidmin, blks.up)
            blks.down = np.where(blks.fodown == blks.fomax, blks.iidmax, blks.down)
            blks.is_fragment = blks.is_fragment | (blks.up > blks.iidmin) | (blks.down < blks.iidmax)

        # If type is to be ignored, all analysis should focus on internal_ids
        elif (fttype == 'any'):
            print(f'Unsuported fttype {fttype}', file=sys.stderr)
            return pd.DataFrame()

        # Unkown fttype
        else:
            print(f'Unknown fttype {fttype}', file=sys.stderr)
            return pd.DataFrame()

        # Return a summary of all regions
        return blks

    def neighbors(self, targets=['query == 1'], before=3, after=3, min_block_distance=0, strand=None, fttype='same', min_block_id=0):
        """
        Find sets of rows, representing genomic regions, that are located near a set of targets.

        The user may choose a set of rows (targets) as anchors, whose neighbors will be
        evaluated by user-defined parameters, such as feature type, strand or distance.

        The returned value is copy of a slice of the original NeighborhoodDF
        object with updated block_id, origin, is_fragment and query columns.

        Arguments:
         - targets : (list of) boolean pd.Series or rules to select targets.

            Rules will used identify target rows using Pandas's eval method.

            If all rules return True for a given row, the row is a target and
            its neighbors will be identified.

         - before : keep at most this number of features, of the same type as the target,
                    before each target
         - after  : keep at most this number of features, of the same type as the target,
                    after each target

         - min_block_distance : minimum distance between two consecutive blocks
                                Blocks separated by more features 

         - strand : how to evaluate rows concerning the value of the strand column
                    Possible values for this option are:
                     - None : ignore strand
                     - same : same strand as the targets
                     -    + : positive strand features and targets only
                     -    - : negative strand features and targets only

         - fttype : how to process feature types of neighbors
                    Supported values:
                     - same : consider only features of the same type as the target
                     - any  : ignore feature type and count all features when
                              setting neighborhood boundaries

         - min_block_id : starting number for block_ids, useful if calling this method
                          eqeuntially through multiple rotifer.genome.data.NeighborhoodDF
                          dataframes
        """
        import sys
        import numpy as np

        # Build boolean pandas.Series to mark targets
        select = True
        if not isinstance(targets, list):
            targets = [targets]
        for code in targets:
            if isinstance(code,str):
                try:
                    select &= self.eval(code)
                except:
                    print(f'Rule {code} failed', file=sys.stderr)
            elif isinstance(code,pd.Series):
                select &= code

        # Find targets and their neighborhood
        if (not isinstance(select,pd.Series)) or (select.sum() == 0):
            print(f'No anchors to search for neighbors were found! Revise your list of targets!')
            return pd.DataFrame()

        # Initialize dataframe for each region (blocks)
        select = self.filter(['assembly','nucleotide','internal_id']).assign(query=select)
        blks = self.vicinity(select['query'], before, after, min_block_distance, fttype, min_block_id)
        if blks.empty:
            return pd.DataFrame()

        # Expand list of features within each neighborhood, tag queries
        # IMPORTANT: row number may increase if min_block_distance < 0!
        blks['internal_id'] = pd.Series(blks[['up','down']].values.tolist()).apply(lambda x: range(x[0],x[1]+1))
        blks = blks.filter(['assembly','nucleotide','internal_id','block_id','rid','origin','is_fragment']).explode('internal_id')
        blks = blks.merge(select, on=['assembly','nucleotide','internal_id'], how="left")

        # Replace columns in self with those from blks
        cols = self.columns
        copy = self.drop(['query','block_id','origin','is_fragment'], axis=1).sort_values(['assembly','nucleotide','start','end'])
        copy = copy.merge(blks, left_on=['assembly','nucleotide','internal_id'], right_on=['assembly','nucleotide','internal_id'], how='inner')
        copy.sort_values(['assembly','nucleotide','block_id','rid','internal_id'], ascending=[True,True,True,True,True], inplace=True)
        copy['query'] = np.where((~copy['query'].isna()) & copy['query'],1,0).astype(int)
        copy = copy[cols]

        # Evaluate other criteria for defining blocks
        if strand == 'same':
            tmp  = (copy.assembly   == copy.shift(1).assembly)
            tmp &= (copy.nucleotide == copy.shift(1).nucleotide)
            tmp &= (copy.block_id   == copy.shift(1).block_id)
            tmp &= (copy.strand     == copy.shift(1).strand)
            copy['block_id'] = (~tmp).cumsum()

            # Cleanup
            valid = copy.groupby('block_id').agg({'query':'sum'}).query('query > 0').block_id
            copy = copy.query('block_id in @valid')
            copy['block_id'] = (~(copy.block_id == copy.shift(1).block_id)).cumsum()

        # Set within block row number
        copy.reset_index(inplace=True, drop=True)
        copy['internal_id'] = copy.index.to_series()

        # Return
        return copy
