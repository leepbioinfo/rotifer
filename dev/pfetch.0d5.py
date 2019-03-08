#!/usr/bin/env python3

import os
import argparse
import sys
sys.path.insert(0, '/home/kaihami/mymodules')

import threading
from threading import Thread
import rotifer.core.cli as corecli
import rotifer.core.functions as cf
import time
from tqdm import trange
import tqdm
from tempfile import mkstemp
import subprocess
from os.path import expanduser
import yaml
from Bio import Entrez
from multiprocessing import Process
import argcomplete

__version__ = '0.5'
__authors__ = 'Gilberto Kaihami; Aureliano Guedes'

################
### Function ###
################

### Arguments
# TODO change for a better helper
def parse_cli():
    parser = corecli.parser(description ='Parallel esl-sfetch, with epost and efetch')

    parser.add(rconfig=':cli.acc')

    parser.add(long_arg = '--database',
               short_arg = '-db',
               default = 'nr',
               dest = 'database',
               helper = 'Protein fasta DB (nr,call, bacteria, archaea) (Default: nr)'
               )

    parser.add(long_arg ='--api_key',
                        helper = 'Insert a valid NCBI API key or a valid config file username',
                        default = '')

    parser.add(long_arg = '--show',
    helper = 'Show databases',
    action = 'store_true')

    parser.add(
    short_arg ="-v",
    long_arg ='--verbose',
    action = 'store_true',
    helper = "Verbose"
    )

    parser.add(
    short_arg = "-p",
    long_arg = '--progress',
    action = 'store_true',
    helper = "Progress bar",
    )

    parser.add(short_arg = '-t',
               long_arg = '--threads',
                default = 5,
                arg_type = int,
                helper = 'Number pf threads (Default:5)')

    argcomplete.autocomplete(parser)

    args = parser.parse_args()

    return (args)

################
### Checkers ###
################

def check_tmp():
    # create a temporary file
    check = True
    tmp = 'tmp.txt'
    x = 0
    while check:
        if tmp in os.listdir(os.getcwd()):
             x += 1
             tmp = 'tmp'+str(x)+'.txt'
        else:
             check = False
    return tmp

def done(initial, final):
    test = [x for x in open(final).read().split('\n') if ">" in x]
    sys.stderr.write('''
All Fasta sequences were extracted!
Number of accessions:      {0}
Number of final sequences: {1}
Happy? {2}

(\_/)
(-.-)
(> <)

'''.format(initial,
                     len(test),
                     str(initial == len(test))
                    )
         )

def find_entrez(accs, api_key):
    len_size = 500
    sub_accs = [accs[x:x+len_size] for x in range(0,len(accs), len_size)]
    for acc in sub_accs:

        tries = 0
        while tries < 2:
            try:
                to_save = Entrez.efetch(db = 'protein',
                                        rettype = 'fasta',
                                        retmode = 'text',
                                        id = ','.join(acc),
                                        api_key = api_key)
                out = [x for x in to_save.read().splitlines() if x != '']
                print('\n'.join(out))
                sys.stdout.flush()
                time.sleep(1)
                break
            except:
                tries +=1
                time.sleep(1)

def find_sequences(tmp_f, database, output, thread = 5, api_key = '935f49538aaa065207b4582a5cc4fcd59408'):
    os.system("cat %s | parallel -j 18 esl-sfetch %s 2> %s.err {} |tee %s " % (tmp_f, database, output, output))
    to_open = open('{0}.err'.format(output)).read().splitlines()
    missing_accs = [x.split(' ')[1] for x in to_open if x != '']

    n = len(missing_accs) // thread if len(missing_accs) // thread > 0 else 1
    sub_missing_accs = [missing_accs[x:x+n] for x in range(0, len(missing_accs), n)]
    jobs = []

    for x in range(len(sub_missing_accs)):
        p = Process(target = find_entrez, args = (sub_missing_accs[x], api_key ))
        jobs.append(p)
        p.start()

    try:
        for p in jobs:
            p.join()
    except KeyboardInterrupt:
        for p in jobs:
            p.terminate()
        sys.exit(2)

def progress(initial, final):
    try:
        t = 0
        checker = True
        pbar = tqdm.tqdm(total=initial)
        time.sleep(1)
        old = 0
        while checker:
            test = len([x for x in open(final).read().split('\n') if ">" in x])
            if old == test:
                t+=1

            pbar.update(test-pbar.n)
            old = test
            if test == initial:
                checker = False
                break
            if t >=100:
                checker = False
                pbar.close()
                break
            time.sleep(0.1)

        pbar.close()
        return 0
    except:
        pass

def show():
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config')

    try:
        n = open(os.path.join(db_local_path, 'pfetch.config'))
        db_local = yaml.load(n)

    except:
        pass

    db_global_path = '/home/kaihami/mymodules/rotifer/config/pfetch.config'
    db_global = yaml.load(open(db_global_path))

    try:
        print('## Local Database')
        print('## Path: {0}'.format(os.path.join(db_local_path,'pfetch.config')))
        print('Name\tPath')
        for k,v in db_local['database'].items():
            print(k+'\t'+v)
        print()
    except:
        pass

    try:
        print('## Global Database')
        print('## Path: {0}'.format(db_global_path))
        print('Name\tPath')
        for k,v in db_global['database'].items():
            print(k+'\t'+v)
        sys.exit()

    except:
        pass

def loadDB():
    database = cf.loadConfig(':pfetch.database')
    return database

############
### MAIN ###
############

if __name__ == '__main__':
    Entrez.email = 'kaihami@gmail.com'
    # Arguments
    args = parse_cli()

    # Check if API is valid
    # print(args.api_key)
    api_key = ''
    try:
        api_key = cf.loadAPI(args.api_key)

    except: pass

    if api_key == '':
        api_key = args.api_key

    if api_key == '':
        api_key = '935f49538aaa065207b4582a5cc4fcd59408'
        # sys.stderr.write('Insert valid NCBI API Key\n')


    threads = args.threads if args.threads <= 5 else 5
#    checkOut = args.output # mktemp
    # Checking if it is valid
    if args.show:
        show()
        sys.exit(0)
    databases = loadDB()
    if args.database in databases.keys():
        database = databases[args.database]
    else:
        database = args.database
    # Create a tmp file        tmp = args.accession
    acc = args.accession
    fd, path_accession = mkstemp()

    with open(path_accession, 'a') as f:
        f.write('\n'.join(acc))

    fd2, path_output = mkstemp()
#        if args.output == "":
#            args.output = check_tmp()
    # Peform the parallel search, returns a .err file with no hits
    first = Thread(target = find_sequences, args = (path_accession, database, path_output, threads, api_key))
    time.sleep(0.4)

    if args.progress:
        second = Thread(target = progress, args = (len(acc), path_output))
        second.start()

    first.start()
    try:
        first.join()
        if args.progress:
            second.join()

    except KeyboardInterrupt:
        first.terminate()
        if args.progress:
            second.terminate()
        sys.exit(2)

     # pass if error
    if args.verbose:
        done(len(acc), path_output)
