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
    "-o", '--output',
    default= '',
#    help=argparse.SUPPRESS
    help = "output file"
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
                                                  description = 'Get fasta sequence from a list of accession'))
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
    os.system("cat %s | parallel -j 18 esl-sfetch %s 2> %s.err {} |tee %s " % (tmp_f, database,
                                                                        output, output))
    # Search for the missing sequencces using epost/efetch
    os.system("awk '{print $2}' %s |\
            sed '/^\s*$/d' |\
            epost -db protein -format acc 2> /dev/null|\
            efetch -db protein -format fasta 2> /dev/null| tee -a %s " % (args.output + '.err', args.output))

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
    try:
        checkOut = args.output
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
        # Create a tmp file, because idk how to use it in the parallel
        tmp = args.accession
        tmp_f = check_tmp()

        with open(tmp_f, 'a') as f:
            f.write('\n'.join(tmp))

        if args.output == "":
            args.output = check_tmp()

        # Peform the parallel search, returns a .err file with no hits
        first = Thread(target = find_sequences, args = (tmp_f, database, args.output))
        time.sleep(0.4)
        if args.progress:
            second = Thread(target = progress, args = (len(tmp), args.output))
            second.start()
        first.start()
        first.join()
        if args.progress:
            second.join()
         # pass if error
        if args.verbose:
            done(len(tmp), args.output)
        # Remove temporary files
        if checkOut == "":
            os.system('rm '+ args.output+'*')
        if checkOut != "":
            os.system('rm '+ args.output+".err")
        os.remove(tmp_f)
    except KeyboardInterrupt:

        if checkOut == "":
            if os.path.isfile(args.output):
                os.system('rm '+ args.output+'*')
        if checkOut != "":
            if os.path.isfile(args.output):
                os.system('rm '+ args.output+".err")