#!/usr/bin/env python3

# For now, we use this as part of our script's template
import os
import sys
_add_path = [ os.path.join(os.path.dirname(os.path.dirname(__file__)), "lib"),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "lib", "python" + str(sys.version_info.major) + "." + str(sys.version_info.minor), "site-packages")
        ]
for _d in _add_path:
    if os.path.exists(_d):
        sys.path.insert(0,_d)

import threading
from threading import Thread
import rotifer.core.cli as corecli
from Bio import SeqIO
from io import StringIO
import pandas as pd
import numpy as np
__version__ = 0.1
__authors__ = 'Gianlucca Nicastro'

def parse_args():
    parser = corecli.parser( description= \
        'using alignemnt in any format as input and make a sequence consensus')


   # parser.add(':cli.core')
    parser.add(long_arg = '--format',
               short_arg = '-f',
               helper = 'Sequence input format default = fasta',
               dest = 'format',
               action = 'store',
               default = 'fasta',
               arg_type = str)

    parser.add(
               helper = 'input file',
               nargs = '*',
               dest = 'input',
               action = corecli.action.autoload,
               duplicates = False
              )

    parser.add(long_arg = '--threshold',
               short_arg = '-th',
               helper = 'Insert the threshold to clean the aligment sequences',
               dest = 'threshold',
               default = 0.9,
               action = 'store',
               arg_type = float)
    parser.add(long_arg = '--table',
               short_arg = '-t',
               helper = 'Output format as seqrwos format',
               dest = 'table',
               action = 'store_true')

    args = parser.parse_args()
    return(args)


args=parse_args()

infile= '\n'.join(args.input)
seq_type= args.format
threshold = args.threshold

def seq_to_seq_rows(infile, seq_type):
    idd=[]
    seq=[]
    fasta_sequences = SeqIO.parse(StringIO(infile),seq_type)
    for fasta in fasta_sequences:
        idd.append(fasta.id)
        seq.append(str(fasta.seq))
    data={'ID':idd, 'seq':seq}
    z = pd.DataFrame(data=data)
    z['len']=z.seq.str.len()
    return z
def freq_col(df):
    r_matrix = pd.DataFrame(index=['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T','-'])
    for x in range(len(df.T)):
        a = df[x].value_counts().to_frame()
        r_matrix = r_matrix.join(a)
    return(r_matrix/len(df))
def aln_matrix(seq_row_object):
    df = pd.DataFrame()
    for x in range(len(seq_row_object.iloc[0,1])):
        df[x] = seq_row_object.iloc[:,1].str[x]
    return df

def aln_clean(seq_rows_ob,aln_matrix_ob, freq_matrix_ob,threshold ):
    df = freq_matrix_ob.T
    df = aln_matrix_ob[df[df['-'] < threshold].index]
    df['clean']= df[df.columns[0:]].apply(
    lambda x: ''.join(x),
    axis=1
)
    return df.set_index(seq_rows_ob.ID)[['clean']]


a = seq_to_seq_rows(infile, seq_type)
b = aln_matrix(a)
c = freq_col(b)
df = aln_clean(a,b,c, threshold)
if args.table:
        df.to_csv(sys.stdout, sep='\t')
else:
    for x,y in df.iterrows():
        print ('>{}\n{}'.format(y.name, y.clean))
