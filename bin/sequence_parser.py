#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import threading
from threading import Thread
import rotifer.core.cli as corecli
from Bio import SeqIO
import pandas as pd
__version__ = 0.1
__authors__ = 'Gianlucca Nicastro'


def parse_args():
    parser = corecli.parser( description= \
        'Parse sequences or aligments and transform it in \
different formats. Ex: Fasta to seqrow, aligments to ungaped fasta,etc')


   # parser.add(':cli.core')
    parser.add(long_arg = '--format',
               short_arg = '-f',
               helper = 'Sequence input format default = fasta',
               dest = 'format',
               action = 'store',
               default = 'fasta',
               arg_type = str)

    parser.add(long_arg = '--ungaped',
               short_arg = '-u',
               helper = 'Remove gaps from alignments',
               dest = 'ungap',
               action = 'store_true')

    parser.add(long_arg = '--input',
               short_arg = '-i',
               helper = 'input file',
               dest = 'infile',
               action = 'store',
               arg_type = str)

    parser.add(long_arg = '--distribution',
               short_arg = '-d',
               helper = 'Remove gaps from alignments',
               dest = 'hist',
               action = 'store_true')
    parser.add(long_arg = '--bins',
               short_arg = '-b',
               helper = 'Numbers of bins for histogram analysis',
               dest = 'bins',
               default = 10,
               action = 'store',
               arg_type = int)
    parser.add(long_arg = '--table',
               short_arg = '-t',
               helper = 'Output format as seqrwos format',
               dest = 'table',
               action = 'store_true')
    parser.add(long_arg = '--filter',
               short_arg = '-fi',
               helper = 'filter sequence by the limit set by the paramethers left limi, right limit',
               dest = 'filtered',
               action = 'store_true')
    parser.add(long_arg = '--left_limit',
               short_arg = '-ll',
               helper = 'The size of the smaller sequence to be filtered',
               dest = 'left_limit',
               action = 'store',
               default = 0,
               arg_type = int)
    parser.add(long_arg = '--right_limit',
               short_arg = '-rl',
               helper = 'The size of the highest sequence after the filter',
               dest = 'right_limit',
               action = 'store',
               default = sys.maxsize,
               arg_type = int)




    args = parser.parse_args()
    return(args)


args=parse_args()

infile=args.infile
seq_type= args.format
ungaped = args.ungap
hist=args.hist

idd=[]
seq=[]
ung=[]
fasta_sequences = SeqIO.parse(infile,seq_type)
for fasta in fasta_sequences:
    idd.append(fasta.id)
    ung.append(str(fasta.seq.ungap('-')))
    seq.append(str(fasta.seq))
    data={'ID':idd, 'seq':seq, 'ung':ung}
z = pd.DataFrame(data=data)
z['len']=z.ung.str.len()
if args.filtered:
    z=z[(z.len > args.left_limit) & (z.len < args.right_limit)]

if ungaped:
    z = z[['ID', 'ung', 'len']].rename({'ung': 'seq'}, axis=1)
else:
    z = z[['ID', 'seq', 'len']]
z.sort_values('len', ascending=False, inplace=True)
if not hist:
    z = z[['ID', 'seq']]
    if args.table:
        z.to_csv(sys.stdout, sep='\t', index=None)
    else:
        for x,y in z.iterrows():
            print ('>{}\n{}'.format(y.ID, y.seq))
else:
    print('Total proteins: {}'.format(len(z)))
    from ascii_graph import Pyasciigraph
    a = z.len.value_counts().to_frame().reset_index()
    a= a.sort_values('index')
    a['raw_bin'] = pd.cut(a['index'],args.bins,precision=0)
    a['bin'] = a.raw_bin.apply(lambda x : '{} - {}'.format(int(x.left),int(x.right)))
    test = a.groupby('bin').agg({'len':'sum'}).reset_index().apply(tuple, axis=1)
    graph = Pyasciigraph()
    for line in  graph.graph('count \t sequence size', test):
        print(line)
