#!/usr/bin/env python3
import rotifer.core.log as rlog
import pandas as pd
from rotifer.core.functions import loadClasses
import os


# Core class
class database:
    '''
    What define a genome block?
        - nucleotide
        - type
        - start
        - end
        - strand
        - feature_id
        - block_id
    '''

    def __init__(self, sources = [],
                 additional_sources = [],
                 verbose = 0, log_file = ''):
        '''
        '''

        # from rotifer.genome.ch import ClickHouse
        try:
            base = os.path.join(os.path.realpath(os.path.join(os.path.abspath(__file__), '..', '..') ), 'genome')
        except:
            base = os.path.dirname(os.path.realpath(__name__))

        self._verbose = verbose
        self._log_file = log_file
        self.sources =[]
        self.sources.extend(sources)
        self.config = {}


        self._switch_source = loadClasses( base+'.db' )

        if additional_sources:
            for additional_source in additional_sources:
                self._switch_source.update(loadClasses(additional_source))


    def open(self, source, config = {}):
        '''
        source: source type
        config: configuration
        '''
        self.target = self._switch_source[source](config,
                                             verbose = self._verbose,
                                             log_file = self._log_file)
        if source not in self.sources:
            self.sources.append(source)

        self.config.update({source: config})

        # rlog.log({3: f'Source is {source}'},
        #          level = self._verbose,
        #          log_file = self._log_file,
        #          name = __name__)
        #
        # rlog.log({3: f'Loading config is {config}'},
        #          level = self._verbose,
        #          log_file = self._log_file,
        #          name = __name__)

    def submit(self, accs = '',
               input_type =['protein_acc'],
               **parameters):
        '''
        accs:       A list of Accession (check valid input types)
        input_type: Accession input type
                    - Accepted types:
                        - protein_acc
                        - locus
                        - assembly
                        - nucleotide
        block_id: Initial block id (default: 0)
        filterby: Several filters
        above:    How many CDS above (default: maximum)
        below:    How many CDS below (default: maximum)
        '''
        self._accs = accs
        self._input_type = input_type
        self._parameters = parameters
        self.target.submit(accs = self._accs,
                           input_type = input_type,
                           **parameters)
        pass

    def fetch_next(self, how = 'block'):
        # genome, block, nucleotide, assembly
        return self.target.fetch_next(how = how)

    def missing(self):
        return self.target.missing()

    def configuration(self):
        print(self.target.config)

    def submitted_parameters(self, outformat = 'table'):
        print(f'Accs:\t{self._accs[0:10]} ...')
        print(f'Input type\t{" ".join(self._input_type)}')

        try:
            for k,v in self._parameters.items():
                print(f'{k}\t{v}')

        except:
            pass

    def available_sources(self):
        return list(self._switch_source.keys())


if __name__ == '__main__':
    from rotifer.genome.data import NeighborhoodDF
    s = database(sources = ['clickhouse'], verbose = 3)

    config = {'clickhouse': ':genome.clickhouse',
              'gff': {'b':2}
             }

    for source in s.sources:
        s.open(source, config = config[source])

        s.submit(['GCF_000067045.1','GCF_002094795.1','B7H17_RS10940', 'a'], above = 0, below = 0, block_id = 100, input_type = ['assembly', 'nucleotide', 'locus'])
        # s.submit(['WP_127490459.1', 'AAP99025.1', 'AAP98158.1'], input_type = ['protein_accession'])
    # a = next(s.fetch_next('block'))
    # b = NeighborhoodDF(a)
    #
    # d = database('clickhouse', {1:2}, additional_source = '~/mymodules/test/dyn/dyn')
