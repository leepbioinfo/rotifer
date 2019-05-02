#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import threading
from threading import Thread
import rotifer.core.cli as corecli
from Bio import SeqIO
from io import StringIO
import pandas as pd
import numpy as np
__version__ = 0.4
__authors__ = 'Gianlucca Nicastro'
__rdoc__='''
DESCRIPTION:
using alignemnt in any format as input and make a sequence consensus
'''

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
               short_arg = '-t',
               helper = 'Insert list of threshold to make the consensus sequences',
               dest = 'threshold',
               default = [1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6],
               action = 'store',
               nargs='+',
               arg_type = float)

    parser.add(long_arg = '--gap',
               short_arg = '-g',
               helper = 'Show gaps as dots instead of underlines',
               dest = 'gap',
               default =True,
               action = 'store_false')
    args = parser.parse_args()
    return(args)


args=parse_args()

infile= '\n'.join(args.input)
seq_type= args.format
threshold_list = args.threshold
gap = args.gap


aromatic = ['F','Y', 'W', 'H']
alifatic = ['I','V','L']
hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
positive = ['H', 'K', 'R']
negative = [ 'D', 'E']
charged = positive + negative
polar = charged + ['Q', 'N', 'S', 'T','C']
alcohol = ['S','T']
tiny = ['G', 'A', 'S']
small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
all_aa = ['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T']

group_size = pd.DataFrame.from_dict({'a': len(aromatic),
                        'l': len(alifatic),
                        'h': len(hydrophobic),
                        '+': len(positive),
                        '-': len(negative),
                        'c': len(charged),
                        'p': len(polar),
                        'o': len(alcohol),
                        'u': len(tiny),
                        's': len(small),
                        'b': len(big),
                        '.': len(all_aa)
                       }, orient='index').rename({0:'size'}, axis=1)

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
def freq_col(df, gap=gap):
    r_matrix = pd.DataFrame(index=['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T','_'])
    g_matrix=pd.DataFrame(index=['a', 'l', 'h', '+', '-', 'c', 'p', 'o', 'u', 's', 'b', '.'])
    for x in range(len(df.T)):
        a = df[x].value_counts().to_frame()
        if gap:
            r_matrix = r_matrix.join(a.rename(index={'-': '_'}))
        else:
            r_matrix = r_matrix.join(a.rename(index={'-': '.'}))
        d = {'a':np.where(df[x].isin(aromatic),1,0).sum(),
         'l':np.where(df[x].isin(alifatic),1,0).sum(),
         'h':np.where(df[x].isin(hydrophobic),1,0).sum(),
         '+':np.where(df[x].isin(positive),1,0).sum(),
         '-':np.where(df[x].isin(negative),1,0).sum(),
         'c':np.where(df[x].isin(charged),1,0).sum(),
         'p':np.where(df[x].isin(polar),1,0).sum(),
         'o':np.where(df[x].isin(alcohol),1,0).sum(),
         'u':np.where(df[x].isin(tiny),1,0).sum(),
         's':np.where(df[x].isin(small),1,0).sum(),
         'b':np.where(df[x].isin(big),1,0).sum(),
         '.':df[x].count()
        }
        g_matrix = g_matrix.join(pd.DataFrame.from_dict(d,orient='index').rename({0:x}, axis=1))
    return(r_matrix/len(df), g_matrix/len(df))
def consensus(a,c, c_freq):
    residue=[]
    freq=[]
    group=[]
    g_freq=[]
    for x in range(len(a.T)):
        aa = a[x].to_frame().sort_values(x, ascending=False)
        residue.append(aa.iloc[0].name)
        freq.append(aa.iloc[0][x])

        cc = c[x].to_frame()
        cc=cc[cc[x] >= c_freq].join(group_size).sort_values(x, ascending=False).drop_duplicates(
            subset='size').sort_values('size').iloc[0]
        group.append(cc.name)
        g_freq.append(cc.iloc[0])
    b = pd.DataFrame({'residue':residue, 'frequency': freq, 'group': group, 'g_freq': g_freq})
    b['consensus']= np.where(b.frequency > c_freq, b.residue,np.where(b.g_freq > c_freq,b.group, '.'))
    return "".join(list(b.consensus))

def aln_matrix(seq_row_object):
    df = pd.DataFrame()
    for x in range(len(seq_row_object.iloc[0,1])):
        df[x] = seq_row_object.iloc[:,1].str[x]
    return df

x = seq_to_seq_rows(infile, seq_type)
y = aln_matrix(x)
rm,gm = (freq_col(y))

print('{}\t{}'.format(x.iloc[0][0],x.iloc[0][1]))
for x in threshold_list:
    print ('consensus/{}%\t{}'.format(int(x*100), consensus(rm,gm,x)))
