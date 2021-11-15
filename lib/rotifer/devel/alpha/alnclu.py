#!/usr/bin/env python3

# Input is a dataframe containing 2 hierarchical clusters
# Output is directory named clusters with one alignment per higher cluster per file

import os
import pandas as pd



def alnclu(info, c80e3="c80e3", i1c1="i1c1", index, i=3):
    '''
    Function to generate alignments for all clusters indicated by c80e3
    takes a indexed multifasta by esl-sfetch as index, and aggregates the
    clusters by a higher level cluster, in this case, i1c1.
    parameter i is 3 by default, increasing it exclude small clusters
    index.fa file should be at current directory,
    or it's relative path defined as well as the index.ssi.
    Usage:
    alnclu(pd.Dataframe, 'c8e3', 'i1c1', 'complete.faa')
    '''

curr_dir = os.getcwd()
os.makedirs(f'{curr_dir}/clusters/hhmdb/aln')

for x in info.c80e3.unique():
    if info[info.c80e3 == x].i1c1.nunique() >= i:
        os.mkdir(f'{curr_dir}/clusters/{x}')
        os.chdir(f'{curr_dir}/clusters/{x}')
        info[info.c80e3 == x].i1c1.drop_duplicates().to_csv(f'{x}.{c80e3}.acc', sep="\t", index=None, header=None)
        os.system(f'esl-sfetch -f {curr_dir}/{index} {x}.{c80e3}.acc | mafft --thread 16 --localpair --maxiterate 1000  > {x}.{c80e3}.aln')
    os.chdir(f'{curr_dir}')

os.chdir('clusters')
os.system('for x in */*.aln; do y=$(echo $x|cut -f1 -d"/"); echo \#$y|cat - $x > ./hhmdb/aln/$y.aln;done')

os.chdir('hhmdb')
os.system(f'/home/leep/ggnicastro/bin/build/build_hhdb.sh ./aln ../{index}')

os.chdir('aln')
os.system(f'for x in *.aln; do hhsearch -i $x -d ../{index} -M 50;done')
os.system('python3 /home/leep/ggnicastro/bin/hhsearch_table.py')
