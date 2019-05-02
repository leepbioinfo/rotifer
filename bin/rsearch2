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
    parser.add_argument('-d', '--database',
                        help = 'Protein fasta DB (Default: /databases/fadb/nr.fa)',
                        default='/databases/fadb/nr.fa')

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

    parser.add_argument('--rneighbors',
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

    parser.add_argument('--fun',
                        help = 'Fun',
                        action = 'store_true')

    args = corecli.parseargs(parents = [parser, parser2])

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
            f1.write('\n'.join(acc))
            f1.write('\n')
        return True

    else:
        with open(output+'.acc', 'a') as f:
            f.write('\n'.join(set(input_file)))
            f.write('\n')
        return False

def pfetch(acc, database, output='rsearch'):
    os.system("cat %s | parallel -j 18 esl-sfetch %s {} > %s 2> %s.err" % (output+'.acc', database,
                                output + '.fa', output))

        # Search for the missing sequencces using epost/efetch
    os.system('''awk '{print $2}' %s | sed '/^\s*$/d' |\
              epost -db protein -format acc |\
              efetch -db protein -format fasta >> %s''' % (output+'.err', output+'.fa'))

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
    os.system('''cat {0} |\
            parallel -n 100 -j 10 --pipe --recstart ">" \
            hmmscan {1} - > {2}'''.format(query, db, outname))

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
            sys.stderr.write('Running rpsblast\n')
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
        os.system('~acpguedes/bin/rneighbors -a {0} -b {1} -asm {2} -of table {3} > {4}'.format(
                                    above,
                                    below,
                                    asm,
                                    input_file,
                                    output+'.neighbor'))

    if gacc:
        os.system('~acpguedes/bin/rneighbors -a {0} -b {1} -gacc {2} -of table {3} > {4}'.format(
                                    above,
                                    below,
                                    gacc,
                                    input_file,
                                    output+'.neighbor'))

    else:
        os.system('~acpguedes/bin/rneighbors -a {0} -b {1} -of table {2} > {3}'.format(
                                    above,
                                    below,
                                    input_file,
                                    output+'.neighbor'))

def acc_from_neighbors(ipt, ipt_acc, output):
    os.system('cut -f 7 {0} | tail -n +2 | grep -v {1} | sort -u > {2}'.format(
                                                ipt, ipt_acc, output+'.neighbor.acc'))

def check_pfetch(acc_file, fasta_file):
    acc = open(acc_file, 'r').read().splitlines()
    fasta = [x for x in open(fasta_file, 'r').read().splitlines() if ">" in x]
    sys.stderr.write(str(coretime.time_now()) +' Number of input: ')
    sys.stderr.write(str(len(acc))+ '\n')
    sys.stderr.write(str(coretime.time_now()) +' Number of output: ')
    sys.stderr.write(str(len(fasta)) + '\n')
    if len(acc) == len(fasta):
        sys.stderr.write(str(coretime.time_now())+ ' All fasta were collected from input\n')

    if len(acc) != len(fasta):
        sys.stderr.write()

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
            pfetch(args.accession, args.database, output)
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
                          asm = args.assembly,
                          input_file = output+'.acc',
                          output = output)

            _e = time.time()
            if verbose:
                sys.stderr.write(str(coretime.time_now()) + ' Elapsed time: ')
                sys.stderr.write(coretime.pretty_time(_e - _s) + '\n')

            acc_from_neighbors(output+'.neighbor', output+'.acc',
                               output)
            if verbose:

                sys.stderr.write(str(coretime.time_now()) +' Collecting neighbors fasta\n')
            _s = time.time()
            pfetch(output+'.neighbor.acc', args.database, output+'.neighbor')
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

        #### Join table ####
        if args.rneighbors:
            df = pd.read_csv(output+'.neighbor', sep = '\t')
            results = results1 + results2 # name of resulting tables of archs
            for result in results:
                _ = open(result).read().splitlines()[1:]
                s = {}
                for line in _:
                    x = line.split('\t')
                    s[x[0]] =  x[1]
                    name = result.split('.')[-2]

                df[name] = df['pid'].map(s)
            df.to_csv(output+'.neighbor.tsv', sep = '\t', index = False)




        end = time.time()
        if verbose:
            sys.stderr.write('## Total time: ')
            sys.stderr.write(coretime.pretty_time(end-start)+'\n')

        if args.fun:
            fun()
    except:
        try:
            if verbose:
                sys.stderr.write(coretime.time_now() +
                    ' An error occurred\nCleaning up the mess\n')
                if args.fun:
                    not_fun()
            shutil.rmtree(os.path.join(os.getcwd(), result_dir))

            sys.exit(1)
        except:
            pass
