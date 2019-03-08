#!/usr/bin/env python3

__version__ = '0.01'
__author__ = "Gilberto Kaihami"
import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))

### Import core cli
import rotifer.core.time as coretime
import rotifer.core.cli as corecli
import rotifer.db.neighbors as neighbors
import rotifer.db.cli as sqlcli

# Other imports
import argparse
import pandas as pd
import time
import subprocess
import shutil
from multiprocessing import Pool, Process, Manager

def parse_cli():
    # OK
    parser = corecli.parser()
    parser.add_argument('-db', '--database',
                        help = 'Protein fasta DB (nr, all, bacteria, archaea) (Default: all)',
                        default='all')

    parser.add_argument('accession',
                        help = 'Input file (accession or a fasta file)',
                        nargs = '*',
                        action = corecli.action.autoload, duplicates = False
                        )

    parser.add_argument('--hmmscan',
                        nargs = 1,
                        default=[],
                        help = '''A dictionary containing the output name and the database separated by (=)\n
                        NAME=PATH (Default: pfam=/database/pfam/Pfam)''',
                        action = corecli.action.dict_options)

    parser.add_argument('--rpsblast',
                        nargs = 1,
                        default=[],
                        help = '''A dictionary containing the output name and the database separated by (=)\n
                        NAME=PATH (Default: rps=~rfsouza/data/rpsdb/allprofiles)''',
                        action = corecli.action.dict_options)
    parser.add_argument('--rcluster',
                        default = False,
                        action = 'store_true',
                        help = 'Choose if run rcluster (aka: multimmseqs)')

    parser.add_argument('--rneighbors',
                        default = False,
                        action = 'store_true',
                        help = 'Choose if run rNeighbors')

    parser.add_argument('--run',
                        nargs = 1,
                        help = 'Select hmmscan and/or rpsblast (default: both)',
                        action = corecli.action.extend)

    parser.add_argument('-o', '--output',
                        help = 'Output file name (default: rsearch)',
                        default = 'rsearch')

    parser.add_argument('-v', '--verbose',
                        action = 'store_true',
                        help = 'Show messages')

    parser2 = sqlcli.rneighbors().input(exclude = ['outformat', 'file',
                                                   'debug', 'log'],
                                        add_help = False, no_open = True)

    parser3 = corecli.config().input(exclude = 'fun')

    parser.add_argument('--fun',
                        help = 'Fun',
                        action = 'store_true')

    args = corecli.parseargs(parents = [parser, parser2,parser3], exclude_from_dump = 'accession')

    return args
class rsearch:
    def __init__(self,args):
        self.input_file = args.accession
        self.output = args.output
        self.verbose = args.verbose
        self.run = args.run
        self.hmmscan = args.hmmscan
        self.rpsblast = args.rpsblast
        self.rneighbors = args.rneighbors
        self.above = args.above
        self.below = args.below
        self.genomic_accession = args.genomic_accession
        databases = {'nr': '/databases/fadb/nr.fa',
                    'all': '/databases/fadb/freeze/all.fa',
                    'bacteria': '/databases/fadb/freeze/bacteria.fa',
                    'archaea': '/databases/fadb/freeze/archaea.fa'
                     }

#        /databases/fadb/nr.fa
        if args.database in databases.keys():
            self.database = databases[args.database]
        else:
            self.database = args.database

        assemblies = {'all': '/home/kaihami/genomes/list_genomes/all.ls',
                'bacteria': '/home/kaihami/genomes/list_genomes/bac.ls' ,
                'archaea': '/home/kaihami/genomes/list_genomes/archaea.ls'
                }
        if args.database in assemblies.keys():
            self.assembly = assemblies[args.database]
        else:
            self.assembly = args.assembly

    def makedirs(self, base = 'rsearch'):
        # OK
        '''
        Create a base output directory
        '''
        default = 'rsearch'
        lsdir =  os.listdir()
        x = 1
        while True:
            if default in os.listdir():
                default = 'rsearch'+'_'+str(x)
                x +=1
            try:
                os.mkdir(default)
                break
            except:
                pass
        return os.path.join(os.getcwd(), default)

    def check_fasta(self):
        # OK
        '''
        Check if the input file is a fasta file or accession list
        '''
        s = ''.join(self.input_file[0:1000])
        if '>' in s:
            acc = [x.split(' ')[0].replace('>','') for x in input_file if '>' in x]
            with open(self.output+'.fa', 'a') as f:
                f.write('\n'.join(self.input_file))
                f.write('\n')
            with open(self.output+'.acc', 'a') as f1:
                f1.write('\n'.join(list(set(acc))))
                f1.write('\n')
            return True

        else:
            with open(self.output+'.acc', 'a') as f:
                f.write('\n'.join(set(self.input_file)))
                f.write('\n')
            return False

    def pfetch(self,output = 'rsearch'):
        # OK
        '''
        Get protein fasta from an acc list
        '''
        os.system("cat %s | sort -u |parallel -j 18 esl-sfetch %s {} > %s 2> %s.err" % (self.output+'.acc', self.database,
                                    self.output + '.fa', self.output))

        # Search for the missing sequencces using epost/efetch
        acc_f = open(self.output+'.acc', 'r').read().splitlines()
        fasta = [x.split(' ')[0].replace('>','') for x in open(self.output+'.fa', 'r').read().splitlines() if ">" in x]
        if len(acc_f) != len(fasta):
            os.system('''awk '{print $2}' %s | sed '/^\s*$/d' |\
                      epost -db protein -format acc |\
                      efetch -db protein -format fasta >> %s''' % (self.output+'.err', self.output+'.fa'))
            fasta2 = [x.split(' ')[0].replace('>','') for x in open(self.output+'.fa', 'r').read().splitlines() if ">" in x]
            if len(acc_f) != len(fasta2):
                missings = [x for x in acc_f if x not in fasta2]
                for missing in missings:
                    os.system('efetch -db protein -format fasta -id {0} >> {1} 2> /dev/null'.format(missing, self.output+'.fa'))

    def check_run(self):
        # OK
        '''
        run hmmscan and/or rpsblast
        '''
        if self.run:
            return self.run
        if not self.run:
            return ['hmmscan', 'rpsblast']

    def check_hmmscan(self):
        '''
        check hmmscan
        '''
        if not self.hmmscan:
            return [{'pfam': '/databases/pfam/Pfam'}]
        if hmmscan:
            return self.hmmscan

    def check_rpsblast(self):
        if not self.rpsblast:
            return [{'aravind': '~rfsouza/data/rpsdb/allprofiles'}]
        if self.rpsblast:
            return self.rpsblast

    def rneighborsRun(self, above = 3, below = 3, asm = '', gacc = '',
                      input_file = '', output = ''):
        if asm:
            os.system('~acpguedes/bin/rneighbors -H 127.0.0.1 -u rotifer -d genomes -p 5432 -a {0} -b {1} -asm {2} -of table {3} > {4}'.format(
                                        self.above,
                                        self.below,
                                        self.assembly,
                                        self.output+'.acc',
                                        self.output+'.neighbor'))

        if gacc:
            os.system('~acpguedes/bin/rneighbors -H 127.0.0.1 -u rotifer -d genomes -p 5432 -a {0} -b {1} -gacc {2} -of table {3} > {4}'.format(
                                        self.above,
                                        self.below,
                                        self.genomic_accession,
                                        self.output+'.acc',
                                        self.output+'.neighbor'))

        else:
            os.system('~acpguedes/bin/rneighbors -H 127.0.0.1 -u rotifer -d genomes -p 5432 -a {0} -b {1} -of table {2} > {3}'.format(
                                        self.above,
                                        self.below,
                                        self.output+'.acc',
                                        self.output+'.neighbor'))

    def acc_from_neighbors(self, df2, fasta_file):
        # Get neighbors acc
        fasta = [x.split(' ')[0].replace('>', '') for x in open(fasta_file, 'r').read().splitlines() if ">" in x]
        accs = list(set(df2[~df2['pid'].isin(fasta)]['pid'].values))
        with open(self.output+'.neighbor.acc', 'w') as f:
            f.write('\n'.join(accs))

    def check_pfetch(self, acc_file, fasta_file):
        # Check pfetch result
        acc = open(acc_file, 'r').read().splitlines()
        fasta = [x for x in open(fasta_file, 'r').read().splitlines() if ">" in x]
        sys.stderr.write(str(coretime.time_now()) +' Number of input: ')
        sys.stderr.write(str(len(acc))+ '\n')
        sys.stderr.write(str(coretime.time_now()) +' Number of output: ')
        sys.stderr.write(str(len(fasta)) + '\n')
        if len(acc) == len(fasta):
            sys.stderr.write(str(coretime.time_now()) + ' All fasta were collected from input\n')

        if len(acc) != len(fasta):
            sys.stderr.write(str(coretime.time_now()) + 'The input is different from output number\n')

    def check_rneighbors(acc_file, rneighbors_file):
        pass
    def search(self):
        start = time.time()
        result_dir = self.makedirs()
        os.chdir(os.path.join(os.getcwd(), result_dir))
        run = self.check_run() # return list to be used
        result1 = []
        result2 = []
        if not self.check_fasta():
            if self.verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Collecting fasta\n')
                _s = time.time()
                # continue as Accession
                number = 0
                self.pfetch(self.output)
                if self.verbose:
                    self.check_pfetch(self.output+'.acc', self.output+'.fa')

                _e = time.time()
                if self.verbose:
                    sys.stderr.write(str(coretime.time_now()) + ' Elapsed time: ')
                    sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')
        tables = []
        hmmscan = self.check_hmmscan()
        rpsblast = self.check_rpsblast()

        ###
        jobs = []
        manager = Manager()

        ##
        tables_hmmscan = manager.list()
        tables_rpsblast = manager.list()
        for r in run:
            if 'hmmscan' in run:
                if self.verbose:
                    sys.stderr.write(str(coretime.time_now())+ ' Running hmmscan\n')
                for x in range(len(hmmscan)):
                    for prefix, db in hmmscan[x].items():
                        if self.verbose:
                            sys.stderr.write(str(prefix) + ' : ' + str(db) +'\n\n')
                        p = Process(target = self.hmmscanRun, args = (self.output, tables_hmmscan, prefix, db))
                        jobs.append(p)
                        p.start()
            if 'rpsblast' in run:
                if self.verbose:
                    sys.stderr.write(' Running rpsblast\n')
                for y in range(len(rpsblast)):
                    for prefix, db in rpsblast[y].items():
                        print(db)
                        if self.verbose:
                            sys.stderr.write(str(prefix) + ' : ' + str(db) +'\n\n')
                        p = Process(target = self.rpsblastRun, args = (self.output, tables_rpsblast, prefix, db))
                        jobs.append(p)
                        p.start()

        #ADD RNEIGHBORS
        if self.rneighbors:
            _s = time.time()
            if self.verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Running rneighbors\n')
            _s = time.time()
            p = Process(target = self.rneighborsRun)
            jobs.append(p)
            p.start()
            _e = time.time()
            if self.verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Elapsed time: ')
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')

        #ADD MMSEQS
        if args.rcluster:
            _s = time.time()
            if self.verbose:
                sys.stderr.write(str(coretime.time_now() + ' Starting rcluster\n'))
            p = Process(target = self.rcluster_query)
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()

        if self.rneighbors:
            df = pd.read_csv(self.output+'.neighbor', sep = '\t') # Load rneighbors table outformat
            df2 = df[df['type'] == 'CDS']
            # TODO Check
            self.acc_from_neighbors(df2, self.output+'.fa')
            if self.verbose:
                sys.stderr.write(str(coretime.time_now()) +' Collecting neighbors fasta\n')
            _s = time.time()
            # TODO check
            self.pfetch(self.output+'.neighbor')

            if self.verbose:
                self.check_pfetch(self.output+'.neighbor.acc' , self.output+'.neighbor.fa')
            jobs2 = []


            if self.verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Collecting neighborhood models\n')
            _e = time.time()
            for r in run:
                if 'hmmscan' in run:
                    if self.verbose:
                        sys.stderr.write(str(coretime.time_now())+ ' Running hmmscan\n')
                    for x in range(len(hmmscan)):
                        for prefix, db in hmmscan[x].items():
                            if self.verbose:
                                sys.stderr.write(str(prefix) + ' : ' + str(db) +'\n\n')
                            p = Process(target = self.hmmscanRun, args = (self.output+'.neighbor', tables_hmmscan, prefix, db))
                            jobs2.append(p)
                            p.start()
                if 'rpsblast' in run:
                    if self.verbose:
                        sys.stderr.write(' Running rpsblast\n')
                    for y in range(len(rpsblast)):
                        for prefix, db in rpsblast[y].items():
                            print(db)
                            if self.verbose:
                                sys.stderr.write(str(prefix) + ' : ' + str(db) +'\n\n')
                            p = Process(target = self.rpsblastRun, args = (self.output+'.neighbor', tables_rpsblast, prefix, db))
                            jobs2.append(p)
                            p.start()
            if self.verbose:
                sys.stderr.write(str(coretime.time_now()) +' Elapsed time: ')
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')
            # neighbors
            if args.rcluster:
                p = Process(target = self.rcluster_neighbor)
                jobs2.append(p)
                p.start()
                _e = time.time()
                if self.verbose:
                    sys.stderr.write(str(coretime.time_now() + ' Elapsed time: '))
                    sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')
            for p in jobs2:
                p.join()

            _hmmscan_res = {}
            for res in tables_hmmscan:
                name = res.split('.')[-2]
                if name not in _hmmscan_res.keys():
                    _hmmscan_res[name] = {}
                for line in open(res).read().splitlines()[1:]:
                    x = line.split('\t')
                    _hmmscan_res[name][x[0]] = x[1]
            for db in _hmmscan_res.keys():
                df[db] = df['pid'].map(_hmmscan_res[db])

            _rpsblast_res = {}
            for res in tables_rpsblast:
                name = res.split('.')[-2]
                if name not in _rpsblast_res.keys():
                    _rpsblast_res[name] = {}
                for line in open(res).read().splitlines()[1:]:
                    x = line.split('\t')
                    _rpsblast_res[name][x[0]] = x[1]
            for db in _rpsblast_res.keys():
                df[db] = df['pid'].map(_rpsblast_res[db])

            cols = df.columns.tolist()
            cols = cols[0:7] + cols[16:]+ cols[7:16]
            df = df[cols]
            if args.rcluster:
                dfQuery = pd.read_csv('query.rawtable', sep = '\t')
                dfQuery = dfQuery.drop(columns=['Hierarchical cluster'], axis = 1)
                query_columns = ['pid'] + ['Query.'+str(x) for x in list(dfQuery.columns[1:])]
                dfQuery.columns = query_columns
                df = pd.merge(df, dfQuery, how='left', on = 'pid')
                try:
                    dfNeighbor = pd.read_csv('neighbors.rawtable', sep = '\t')
                    dfNeighbor = dfNeighbor.drop(columns=['Hierarchical cluster'], axis = 1)
                    neighbor_columns = ['pid'] + ['Neighbors.'+str(x) for x in list(dfNeighbor.columns[1:])]
                    dfNeighbor.columns = neighbor_columns
                    df = pd.merge(df, dfNeighbor, how='left', on = 'pid')

                    filtered_df = df[df[query_columns[1]].notnull()][query_columns[0:2]]
                    query_filter = dict(zip(filtered_df['pid'], filtered_df[query_columns[1]].astype(str)+'.Query' ) )
                    filtered_df2 = df[df[neighbor_columns[1]].notnull()][neighbor_columns[0:2]]
                    neighbor_filter = dict(zip(filtered_df2['pid'], filtered_df2[neighbor_columns[1]].astype(str)+'.Neighbor' ) )
                    query_filter.update(neighbor_filter)

                    df['Cluster'] = df['pid'].map(query_filter)
                except:
                    pass

            df.to_csv(self.output+'.neighbor.tsv', sep = '\t', index = False)


    def rcluster_query(self):
        os.system('~kaihami/bin/multimmseqs -i {0} -c 0.8 --identity 0.3 0.4 0.5 0.6 -b query > /dev/null 2> /dev/null'.format(self.output+'.fa'))
    def rcluster_neighbor(self):
        os.system('~kaihami/bin/multimmseqs -i {0} -c 0.6 --identity 0.3 0.4 0.5 0.6 -b neighbors > /dev/null 2> /dev/null'.format(self.output+'.neighbor.fa'))

    def hmmscanRun(self, query, tables, prefix, db):
        outname = query+'.'+prefix
        ### Use process call error
        os.system('hmmscan --cpu 10 {0} {1} > {2}'.format(db, query + '.fa', outname))
        os.system('hmmer2table -s {0} | domain2architecture > {1}'.format(query+'.'+prefix,
                                                                          query+'.'+prefix+'.arch'))
        return tables.append(query+'.'+prefix+'.arch')

    def rpsblastRun(self, query, tables, prefix, db):
        outname = query+'.'+prefix
        ### Use process call error
        os.system('''rpsblast -num_threads 10 \
                    -db {0} -query {1} > {2}
                  '''.format(db, query + '.fa', outname))
        os.system('blast2table -s {0} | domain2architecture > {1}'.format(query+'.'+prefix,
                                                                          query+'.'+prefix+'.arch'))
        return tables.append(query+'.'+prefix+'.arch')

def fun():
    sys.stderr.write('''
    (\_/)
    (^_^)
    (")(") \n''')

def not_fun():
    sys.stderr.write('''
         (\_/)
         (. .)
        C(")(") \n''')


if __name__ == '__main__':
    try:
        args = parse_cli()
        rsearch(args).search()
    except (KeyboardInterrupt):
#        shutil.
        print("Stop me!")
        sys.exit(0)
