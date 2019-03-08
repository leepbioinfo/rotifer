#!/usr/bin/env python3

import os
import argparse
import sys
sys.path.insert(0, '/home/kaihami/mymodules')

import threading
from threading import Thread

import rotifer.core.cli as corecli

import time
from tqdm import trange
import tqdm
from tempfile import mkstemp
import subprocess

__version__ = '0.1'
__authors__ = 'Gilberto Kaihami; Aureliano Guedes'

################
### Function ###
################

#def msg(name = None):
#        return '''
#    '''

### Arguments
# TODO change for a better helper
def parse_cli():
    parser = argparse.ArgumentParser(description = 'Parallel esl-sfetch, with epost and efetch',
                                     formatter_class = argparse.RawTextHelpFormatter) #, usage = msg())
    parser.add_argument(
    "-db", '--database',
    help = "Protein fasta DB (nr, all, bacteria, archaea) (Default: all)",
    default = 'all'
#    help=argparse.SUPPRESS
    )

    parser.add_argument(
    'accession',
    help = "Input acession file",
#    help=argparse.SUPPRESS,
    nargs = '*',
    action = corecli.action.autoload
    )

    parser.add_argument(
    "-v", '--verbose',
    action = 'store_true',
    #    help=argparse.SUPPRESS
    help = "Verbose"
    )
    parser.add_argument(
    "-p", '--progress',
    action = 'store_true',
    #    help=argparse.SUPPRESS
    help = "Progress bar"
    )
    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'pfetch',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Get fasta sequence from a list of accession'
                                                  ))
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
def find_sequences(tmp_f, database, output):
    os.system("cat %s | parallel -j 18 esl-sfetch %s 2> %s.err {} |tee %s " % (tmp_f, database, output, output))
#    process = subprocess.Popen('cat %s | parallel -j 18 esl-sfetch %s 2> %s.err {}' % (tmp_f, database, output), shell = True, stdout = subprocess.PIPE)
#    c = process.communicate()[0].decode('utf-8')[:-1]
#    print(c)
    # Search for the missing sequencces using epost/efetch
    os.system("awk '{print $2}' %s |\
            sed '/^\s*$/d' |\
            epost -db protein -format acc 2> /dev/null|\
            efetch -db protein -format fasta 2> /dev/null| tee -a %s " % (output+ '.err', output))

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

############
### MAIN ###
############

if __name__ == '__main__':
    # Arguments
    args = parse_cli()
#    checkOut = args.output # mktemp
    # Checking if it is valid

    databases = {'nr': '/databases/fadb/nr.fa',
                'all': '/databases/fadb/freeze/all.fa',
                'bacteria': '/databases/fadb/freeze/bacteria.fa',
                'archaea': '/databases/fadb/freeze/archaea.fa'
                 }
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
    first = Thread(target = find_sequences, args = (path_accession, database, path_output))
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
