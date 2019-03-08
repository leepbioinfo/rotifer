#!/usr/bin/env python3
import pandas as pd
import rotifer.core.log as rlog
from rotifer.core.functions import loadConfig

class clickhouse:
    def __init__(self, config= {},
                 table_name = '',
                 uri = '',
                 verbose = 0,
                 log_file =''):

        '''
        Connect to clickhouse
        '''

        from rotifer.alchemy.connect import clickhouse

        rlog.log({1:'click'}, level = 1, name = __name__)
        if config.startswith(':'):
            config = (loadConfig(config))

        self.config = config
        self.uri = 'clickhouse://default:@localhost/rotifer'
        self.table_name = 'genomes2'
        self._log_file = log_file
        self._verbose = verbose

        for c in config.keys():
            if c in ['uri', 'table_name']:
                try:
                    setattr(self,c, config[c])
                except:
                    if self.uri == '' or self.table_name == '':
                        rlog.log({2: f'Argument not {c} is empty'
                                  },
                            level = self._verbose,
                            log_file = self._log_file,
                            name = __name__)

            else:
                rlog.log({2: f'Argument not used in the config[{c}]',
                          3: f'Loading CH configuration argument not used [{c}]'
                          },
                    level = self._verbose,
                    log_file = self._log_file,
                    name = __name__)

        try:
            rlog.log({3: f'Connecting to ClickHouse DB',
                      },
                level = self._verbose,
                log_file = self._log_file,
                name = __name__)

            db =  clickhouse(uri = self.uri,
                                    table_name = self.table_name)
            self._conn = db.conn
            self._db_columns = db.table_columns()
            rlog.log({3: f'Connected to ClickHouse DB',
                      },
                level = self._verbose,
                log_file = self._log_file,
                name = __name__)

        except:
            rlog.log({
                 3: f'Unable to connect to ClickHouse Database'
                },
                level = self._verbose,
                log_file = self._log_file,
                name = __name__
                )

    def submit(self, accs,
               input_type = ['protein_acc'],
               **parameters):

        '''
        Submit/commit a query
        ----------
        accs:       A list of accessions (check valid input types)
        input_type: Accession input type
                        - Accepted types:
                            - protein_acc
                            - locus
                            - assembly
                            - nucleotide
                            - organism
                            - gene

        parameters:
            block_id:
            filterby:
            above:
            below:

        '''

        _switch_ipttype = {
                           'protein_acc': 'pid',
                           'locus': 'locus',
                           'assembly': 'assembly',
                           'nucleotide': 'nucleotide',
                           'organism': 'organism',
                           'gene': 'gene'
                          }

        try:
            self._filterby = parameters['filterby']
        except:
            self._filterby = {}

        try:
            self._above = parameters['above']
        except:
            self._above = 'max'

        try:
            self._below = parameters['below']
        except:
            self._below = 'max'

        try:
            self._distance = parameters['distance']
        except:
            self._distance = 0

        try:
            self._block_id = parameters['block_id']
        except:
            self._block_id = 0

        # Abstraction of the input type
        # Convert an arbitrary input type to a valid ipt type
        self._input_type = []
        for ipt in input_type:
            try:
                self._input_type.append(_switch_ipttype[ipt])

            except:
                rlog.log({
                     2: f'ClickHouse invalid input type [{ipt}]',
                     3: f'Invalid input type for ClickHouse [{ipt}], valid options [{", ".join(list(_switch_ipttype.keys()) )}]'
                    },
                    level = self._verbose,
                    log_file = self._log_file,
                    name = __name__
                    )

        self.accs = accs

        # Testing
        if accs == '':
            self.accs = ['XP_753543.1', 'OXN08776.1']
        else:
            self.accs = accs

        # self._input_type = inputtype

        self._nucleotide = self._find_nucletides()


    def fetch_all(self):
        # Get all hits
        pass

    def fetch_next(self, how, **kwargs):

        # genome, assembly, block
        # Generator
        _switch_fetch = {
                        'block': self._next_block,
                        'nucleotide': self._next_nucleotide,
                        'assembly': self._next_assembly,
                        'genome': self._next_assembly
                         }

        return _switch_fetch[how]()

    def fetch_one(self):
        # Get first hit
        pass

    def configuration(self):
        # Show database configuration
        print(self.config)

    def missing(self):
        res = []
        for ipt in self._input_type:
            for x in range(0, len(self._input_type)):
                q_ipt_type = 'WHERE '
                s = self._conn.execute(f'''SELECT DISTINCT {ipt}
                                            FROM {self.table_name}
                                            WHERE {ipt} IN
                                             ({', '.join(self.accs_formated)})
                                             ''').fetchall()
                if s:
                    _ = list(zip(*s))[0]
                    res.extend(_)
                    break

        # List comprehension find missing values
        miss = [x for x in self.accs if x not in res]

        return miss
    def _find_nucletides(self):
        # Get assemblies here

        ### BLOCK FOR QUERY
        # Build query for type
        # Check if valid

        q_ipt_type = 'WHERE '

        self.accs_formated = ["'" + x + "'" for x in self.accs]

        for ipt_type in self._input_type:
            if ipt_type in self._db_columns:
                q_ipt_type += f"{ipt_type} IN ({','.join(self.accs_formated)}) OR "
            else:
                rlog.log({2: f'No column named [{ipt_type}]'}, level = 2)

        self.q_ipt_type = q_ipt_type.rstrip(' OR ')

        ### END BLOCK FOR QUERY

        self.nucleotide_in_db = self._conn.execute(f'''SELECT DISTINCT nucleotide, nuc_asm
                                            FROM {self.table_name}
                                            {self.q_ipt_type}
                                             ''')

        return self.nucleotide_in_db

    def _next_nucleotide(self, block_id = 0):
        pass

    def _next_assembly(self, block_id = 0):
        pass

    def _next_block(self):
        # This is a list of tuples
        for nucleotide,nuc_asm in self._nucleotide:

            # Filter here?
            accs_indexes = self._conn.execute(f"""SELECT feature_order FROM {self.table_name}
                                         WHERE nuc_asm = '{nuc_asm}'
                                              AND nucleotide = '{nucleotide}'
                                        {self.q_ipt_type.replace('WHERE','AND')}""")


            max_index, min_index = self._conn.execute(f"""SELECT max(feature_order), min(feature_order) from {self.table_name}
                                            where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                            type = 'CDS'""").fetchone()
            topology = self._conn.execute(f"""SELECT distinct(topology) from {self.table_name}
                                        where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                        type = 'CDS'""").fetchone()[0]

            if self._below == 'max':
                self._below = min_index
            if self._above == 'max':
                self._above = max_index

            columns = ['internal_id','nucleotide',
                       'assembly', 'topology', 'start', 'end',
                       'strand','pid', 'type', 'plen', 'locus',
                       'seq_type', 'gene', 'product','organism',
                       'classification']

            accs_indexes = sorted(list(zip(*accs_indexes))[0])

            intervals = []

            for value in accs_indexes:

                find_up = value - self._above
                find_down = value + self._below
                if intervals:
                    if intervals[-1][1] - find_up+1 >= self._distance:
                        if self._above+self._below+1 >= max_index:
                            intervals[len(intervals)-1] = [0, max_index]
                        else:
                            intervals[len(intervals) -1 ] = [intervals[len(intervals)-1][0], find_down]
                    else:
                        if self._above + self._below +1 >= max_index:
                            intervals.append([0, max_index])
                        else:
                            intervals.append([find_up, find_down])

                else:
                    if self._above +self._below +1 >= max_index:
                        intervals.append([0, max_index])
                    else:
                        intervals.append([find_up, find_down])

            for up_index, down_index in intervals:
                sub_df = pd.DataFrame()
                get_more_up = ''
                get_more_down = ''
                if up_index < min_index:
                    if topology != 'linear':
                        get_more_up = max_index+ up_index
                    up_index = min_index

                if down_index > max_index:
                    if topology != 'linear':
                        get_more_down = down_index - max_index
                    down_index = max_index

                if get_more_up:
                    _ = self._conn.execute(f"""
                                        SELECT internal_id from {self.table_name}
                                        WHERE nuc_asm = '{nuc_asm}' and
                                        nucleotide = '{nucleotide}' and
                                        type = 'CDS'
                                        and feature_order in ({max_index}, {get_more_up})
                                          """).fetchall()
                    _ = [x[0] for x in _]

                    if len (_) == 1:
                        max_order = min_order = _[0]
                    else:
                        max_order, min_order = sorted(_)

                    q = self._query_region(nuc_asm = nuc_asm,
                                           nucleotide = nucleotide,
                                           max_order = max_order,
                                           min_order = min_order)

                    sub_df = pd.DataFrame(q,
                                          columns = columns)

                if not sub_df.empty:
                    _ = self._conn.execute(f"""
                                        SELECT internal_id from {self.table_name}
                                        WHERE nuc_asm = '{nuc_asm}' and
                                        nucleotide = '{nucleotide}' and
                                        type = 'CDS'
                                        and feature_order in ({up_index}, {down_index})
                                          """).fetchall()

                    _ = [x[0] for x in _]

                    if len (_) == 1:
                        max_order = min_order = _[0]
                    else:
                        max_order, min_order = sorted(_)

                    q = self._query_region(nuc_asm = nuc_asm,
                                           nucleotide = nucleotide,
                                           max_order = max_order,
                                           min_order = min_order)

                    t = pd.DataFrame(q,
                                          columns = columns)

                    sub_df = pd.concat([sub_df, t])

                else:
                    _ = self._conn.execute(f"""
                                        SELECT internal_id from {self.table_name}
                                        WHERE nuc_asm = '{nuc_asm}' and
                                        nucleotide = '{nucleotide}' and
                                        type = 'CDS'
                                        and feature_order in ({up_index}, {down_index})
                                          """).fetchall()
                    _ = [x[0] for x in _]
                    if len (_) == 1:
                        max_order = min_order = _[0]
                    else:
                        max_order, min_order = sorted(_)

                    q = self._query_region(nuc_asm = nuc_asm,
                                           nucleotide = nucleotide,
                                           max_order = max_order,
                                           min_order = min_order)

                    sub_df = pd.DataFrame(q,
                                          columns = columns)

                if get_more_down:
                    _ = self._conn.execute(f"""
                                        SELECT internal_id from {self.table_name}
                                        WHERE nuc_asm = '{nuc_asm}' and
                                        nucleotide = '{nucleotide}' and
                                        type = 'CDS'
                                        and feature_order in ({get_more_down})
                                          """).fetchall()
                    try:
                        _ = [x[0] for x in _]
                        if len (_) == 1:
                            max_order = min_order = _[0]
                        else:
                            max_order, min_order = sorted(_)

                        q = self._query_region(nuc_asm = nuc_asm,
                                               nucleotide = nucleotide,
                                               max_order = max_order,
                                               min_order = min_order)

                        if not sub_df.empty:
                            t = pd.DataFrame(q,
                                              columns = columns)

                            sub_df = pd.concat([sub_df, t])
                    except:
                        pass


                if not sub_df.empty:
                    self._block_id +=1
                    sub_df['block_id'] = sub_df.shape[0]*[self._block_id]

                    sub_df['query'] = sub_df['pid'].map(lambda x: 1 if x in self.accs else 0)
                    if self._block_id == 1:
                        print_header = True
                    else:
                        print_header = False

                    self.header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene modified product organism classification'.split()

                    try:
                        dc = {k:v for k,v in original_acc}
                        sub_df['modified'] = sub_df['pid'].map(lambda x: dc[x] if x in dc.keys() else '.')

                    except:
                        sub_df['modified'] = '.'

                    sub_df = sub_df.drop_duplicates()

                    yield sub_df[self.header] # this is a object for data

        self._nucleotide

    def _query_region(self, nuc_asm,
                      nucleotide,
                      max_order,
                      min_order):
        q = (self._conn.execute(f"""SELECT
                                      internal_id,nucleotide,
                                               assembly, topology, start, end,
                                               strand,pid, type, plen, locus,
                                               seq_type, gene, product,organism,
                                               taxonomy
                              FROM {self.table_name} where nuc_asm = '{nuc_asm}'
                                                  and nucleotide = '{nucleotide}'
                                                  and internal_id >= {max_order}
                                                  and internal_id <= {min_order}
                              ORDER BY internal_id""").fetchall())
        return q

#gdb.open(source = source, config = entrezopts)
#s.submit(accs, input_type = ['protein accession'], filterby=args.filterby, after=7, before=7)
# rgd = s.fetch_all()
# for rgd in s.fetch_next("genome"): # fetch next rotifer.genome.data (i.e. collection of dataframes)
# for block in rgd.blocks():

