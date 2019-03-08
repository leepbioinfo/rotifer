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

def parse_cli():
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
                        NAME=PATH (Default: pfam=/database/pfam/PFAM)''',
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

def makedir(base = 'rsearch'):
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

def check_fasta(input_file,base, output = 'rsearch'):
    '''
    Check if the input file is a fasta file or accession list
    '''
    s = ''.join(input_file[0:1000])
    if '>' in s:
        acc = [x.split(' ')[0].replace('>','') for x in input_file if '>' in x]
        with open(output+'.fa', 'a') as f:
            f.write('\n'.join(input_file))
            f.write('\n')
        with open(output+'.acc', 'a') as f1:
            f1.write('\n'.join(list(set(acc))))
            f1.write('\n')
        return True

    else:
        with open(output+'.acc', 'a') as f:
            f.write('\n'.join(set(input_file)))
            f.write('\n')
        return False

def pfetch(acc, database, output='rsearch'):
    os.system("cat %s | sort -u |parallel -j 18 esl-sfetch %s {} > %s 2> %s.err" % (output+'.acc', database,
                                output + '.fa', output))

    # Search for the missing sequencces using epost/efetch
    acc_f = open(output+'.acc', 'r').read().splitlines()
    fasta = [x.split(' ')[0].replace('>','') for x in open(output+'.fa', 'r').read().splitlines() if ">" in x]
    if len(acc_f) != len(fasta):
        os.system('''awk '{print $2}' %s | sed '/^\s*$/d' |\
                  epost -db protein -format acc |\
                  efetch -db protein -format fasta >> %s''' % (output+'.err', output+'.fa'))
        fasta2 = [x.split(' ')[0].replace('>','') for x in open(output+'.fa', 'r').read().splitlines() if ">" in x]
        if len(acc_f) != len(fasta2):
            missings = [x for x in acc_f if x not in fasta2]
            for missing in missings:
                os.system('efetch -db protein -format fasta -id {0} >> {1} 2> /dev/null'.format(missing, output+'.fa'))

def check_run(run):
    if run:
        return run
    if not run:
        return ['hmmscan', 'rpsblast']

def check_hmmscan(hmmscan):
    if not hmmscan:
        return [{'pfam': '/databases/pfam/Pfam'}]
    if hmmscan:
        return hmmscan

def check_rpsblast(rpsblast):
    if not rpsblast:
        return [{'rps': '~rfsouza/data/rpsdb/allprofiles'}]
    if rpsblast:
        return rpsblast

def hmmscanRun(query, db, output='rsearch', prefix = 'pfam'):
    outname = output+'.'+prefix
    ### Use process call error
    os.system('hmmscan --cpu 10 {0} {1} > {2}'.format(db, query, outname))

def rpsblastRun(db, query, output='rsearch', prefix = 'rps'):
    outname = output+'.'+prefix
    ### Use process call error
    os.system('''rpsblast -num_threads 10 \
                -db {0} -query {1} > {2}
              '''.format(db, query, outname))

def rsearchRun(run, hmmscan, rpsblast, query, output, verbose):
    tables = []
    if 'hmmscan' in run:
        if verbose:
            sys.stderr.write(str(coretime.time_now())+ ' Running hmmscan\n')
        for scan in hmmscan:
            for prefix, db in scan.items():
                if verbose:
                    sys.stderr.write(str(prefix) + ' : ' + str(db) +'\n\n')
                hmmscanRun(query, db, output, prefix)
                os.system('hmmer2table -s {0} | domain2architecture > {1}'.format(output+'.'+prefix,
                                                                                  output+'.'+prefix+'.arch'))
                tables.append(output+'.'+prefix+'.arch')

    if 'rpsblast' in run:
        if verbose:
            sys.stderr.write(' Running rpsblast\n')
        for scan in rpsblast:
            for prefix, db in scan.items():
                if verbose:
                    sys.stderr.write(str(prefix) + ' : ' + str(db) +'\n\n')
                rpsblastRun(db, query, output, prefix)
                os.system('blast2table -s {0} | domain2architecture > {1}'.format(output+'.'+prefix,
                                                                                  output+'.'+prefix+'.arch'))
                tables.append(output+'.'+prefix+'.arch')

    return tables

def rneighborsRun(above = 3, below = 3, asm = '', gacc = '',
                  input_file = '', output = ''):
    if asm:
        os.system('~acpguedes/bin/rneighbors -H 127.0.0.1 -u rotifer -d genomes -p 5432 -a {0} -b {1} -asm {2} -of table {3} > {4}'.format(
                                    above,
                                    below,
                                    asm,
                                    input_file,
                                    output+'.neighbor'))

    if gacc:
        os.system('~acpguedes/bin/rneighbors -H 127.0.0.1 -u rotifer -d genomes -p 5432 -a {0} -b {1} -gacc {2} -of table {3} > {4}'.format(
                                    above,
                                    below,
                                    gacc,
                                    input_file,
                                    output+'.neighbor'))

    else:
        os.system('~acpguedes/bin/rneighbors -H 127.0.0.1 -u rotifer -d genomes -p 5432 -a {0} -b {1} -of table {2} > {3}'.format(
                                    above,
                                    below,
                                    input_file,
                                    output+'.neighbor'))

def acc_from_neighbors(df2, fasta_file, output):
    fasta = [x.split(' ')[0].replace('>', '') for x in open(fasta_file, 'r').read().splitlines() if ">" in x]
    accs = list(set(df2[~df2['pid'].isin(fasta)]['pid'].values))
    with open(output+'.neighbor.acc', 'w') as f:
        f.write('\n'.join(accs))

def check_pfetch(acc_file, fasta_file):
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
        input_file = args.accession
        output = args.output
        verbose = args.verbose

        databases = {'nr': '/databases/fadb/nr.fa',
                    'all': '/databases/fadb/freeze/all.fa',
                    'bacteria': '/databases/fadb/freeze/bacteria.fa',
                    'archaea': '/databases/fadb/freeze/archaea.fa'
                     }
#        /databases/fadb/nr.fa
        if args.database in databases.keys():
            database = databases[args.database]
        else:
            database = args.database

        assemblies = {'all': '/home/kaihami/genomes/list_genomes/all.ls',
                'bacteria': '/home/kaihami/genomes/list_genomes/bac.ls' ,
                'archaea': '/home/kaihami/genomes/list_genomes/archaea.ls'
                }
        if args.database in assemblies.keys():
            assembly = assemblies[args.database]
        else:
            assembly = args.assembly
        # Time
        start = time.time()

        # Create a result dir
        result_dir = makedir()
        os.chdir(os.path.join(os.getcwd(), result_dir))

        # Decide run => return a list
        run = check_run(args.run)

        # Check if input is a fasta file
        fasta = check_fasta(input_file, result_dir, output = output)

        # Check hmmscan/rpsblast output and database
        hmmscan = check_hmmscan(args.hmmscan)
        rpsblast = check_rpsblast(args.rpsblast)
        result1 = []
        result2 = []
        # Check if the input is a fasta file
        if not fasta:
            if verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Collecting fasta\n')
            _s = time.time()
            # continue as Accession
            number = 0
            pfetch(args.accession, database, output)
            if verbose:
                check_pfetch(output+'.acc', output+'.fa')

            _e = time.time()
            if verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Elapsed time: ')
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')

        query = output+'.fa'
        results1 = rsearchRun(run, hmmscan, rpsblast, query, output, verbose)

        if args.rneighbors:
            if verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Running rneighbors\n')
            _s = time.time()
            rneighborsRun(above= args.above, below = args.below,
                          gacc = args.genomic_accession,
                          asm = assembly,
                          input_file = output+'.acc',
                          output = output)

            _e = time.time()
            if verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Elapsed time: ')
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')

            df = pd.read_csv(output+'.neighbor', sep = '\t') # Load rneighbors table outformat
            df2 = df[df['type'] == 'CDS']
            # Check if right
            acc_from_neighbors(df2, output+'.fa',
                               output)

            if verbose:

                sys.stderr.write(str(coretime.time_now()) +' Collecting neighbors fasta\n')
            _s = time.time()
            pfetch(output+'.neighbor.acc', database, output+'.neighbor')

            if verbose:
                check_pfetch(output+'.neighbor.acc' , output+'.neighbor.fa')

            if verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Collecting neighborhood models\n')
            results2 = rsearchRun(run, hmmscan, rpsblast, output+'.neighbor.fa', output+'.neighbor', verbose)

            _e = time.time()
            if verbose:
                sys.stderr.write(str(coretime.time_now()) +' Elapsed time: ')
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')

        ### MMSeqs here
        if args.rcluster:
            # query
            _s = time.time()
            if verbose:
                sys.stderr.write(str(coretime.time_now() + ' Starting rcluster\n'))
            os.system('~kaihami/bin/multimmseqs -i {0} -c 0.8 --identity 0.3 0.4 0.5 0.6 -b query > /dev/null 2> /dev/null'.format(output+'.fa'))
            # neighbors
            os.system('~kaihami/bin/multimmseqs -i {0} -c 0.6 --identity 0.3 0.4 0.5 0.6 -b neighbors > /dev/null 2> /dev/null'.format(output+'.neighbor.fa'))
            _e = time.time()
            if verbose:
                sys.stderr.write(str(coretime.time_now() + ' Elapsed time: '))
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')

        #### Join table ####
        if args.rneighbors:
            results = results1 + results2 # name of resulting tables of archs
            for y in range(len(results1)):
                res1 = open(results1[y]).read().splitlines()[1:]
                res2 = open(results2[y]).read().splitlines()[1:]
                s = {}
                for line in res1:
                    x = line.split('\t')
                    s[x[0]] =  x[1]
                name = results1[y].split('.')[-2]
                print(results1[y])
                for line in res2:
                    x = line.split('\t')
                    s[x[0]] =  x[1]
                df[name] = df['pid'].map(s)

            # reorder columns

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

            df.to_csv(output+'.neighbor.tsv', sep = '\t', index = False)

        end = time.time()
        if verbose:
            sys.stderr.write('## Total time: ')
            sys.stderr.write(coretime.pretty_time(end-start)+'\n')

        if args.fun:
            fun()
    except (KeyboardInterrupt, SystemExit):
#        shutil.
        print("Stop me!")
        sys.exit(0)
