#!/usr/bin/env python3
import subprocess
import argparse
import sys
import rotifer.core.cli as corecli
import tqdm

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Get fasta sequence')
    parser.add_argument(
    'locus',
    help = "Input locus tag file",
    nargs = '*',
    action = corecli.action.autoload
    )
    parser.add_argument(
    '-p', '--progress',
    action = 'store_true',
    help = 'Progress bar'
    )
    parser.add_argument(
    '-a', '--annotation',
    action = 'store_true',
    help = 'Insert original locus tag using a Header ###')
    parser.add_argument(
    '-t', '--table',
    action = 'store_true',
    help = 'Output table format (Locus, ID, Header, Seq)')
    args = parser.parse_args()
    return args

def fasta2line(d, locus):
    tmp = [x for x in d.split('\n') if x != '']
#    print(tmp)
    dc = {}
    for x in range(0, len(tmp)):
        if tmp[x].startswith('>'):
            header = tmp[x]
            dc[header] = ''
        else:
            dc[header] += tmp[x]
    ls = []
    for k, v in dc.items():
        print('\t'.join([locus, k.split(' ')[0].replace('>',''),
                   k,
                   v]))

if __name__ == '__main__':
    args = parse_cli()
    loci = args.locus
    if args.table:
        header = ['Loucs', 'ID','Header', 'Seq']
        print('\t'.join(header))
    if args.progress:
        pbar = tqdm.tqdm(total = len(loci))
        try:
            for locus in loci:
                if args.annotation:
                    print('###', locus)
                p = subprocess.Popen(['esearch -db protein -query "{0}" | efetch -format fasta'.format(locus)],
                                     shell = True, stdout=subprocess.PIPE)
                com = p.communicate()[0]
                d = com.decode('utf-8')
                if args.table:
                    fasta2line(d, locus)
                else:
                    print(d) # print fasta
                sys.stdout.flush()
                pbar.update(1)
        except KeyboardInterrupt:
            raise Exception('Hey stop!')
        pbar.close()
    else:
        try:
            for locus in loci:
                if args.annotation:
                    print('###', locus)
                p = subprocess.Popen(['esearch -db protein -query "{0}" | efetch -format fasta'.format(locus)],
                                     shell = True, stdout=subprocess.PIPE)
                com = p.communicate()[0]
                d = com.decode('utf-8')
                if args.table:
                    fasta2line(d, locus)
                else:
                    print(d) # print fasta
                sys.stdout.flush()
        except KeyboardInterrupt:
            raise Exception('Hey stop!')

