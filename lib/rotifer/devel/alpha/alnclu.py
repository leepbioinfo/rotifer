#!/usr/bin/env python3

# Input is a dataframe containing 2 hierarchical clusters
# Output is directory named clusters with one alignment per higher cluster per file


def alnclu(info, c80e3="c80e3", c100i100="c80i70", i=3, index="complete.fa"):
    '''
    Function to generate alignments for all clusters indicated by c80e3
    takes a indexed multifasta by esl-sfetch as index, and aggregates the
    clusters by a higher level cluster, in this case, c100i100.
    parameter i is 3 by default, increasing it exclude small clusters
    index.fa file should be at current directory,
    or it's relative path defined as well as the index.ssi.
    Usage:
    alnclu(pd.Dataframe, 'c80e3', 'c100i100', 'complete.faa')
    '''
    import os
    import pandas as pd
    import rotifer.devel.beta.sequence as rdbs

    # Validate input
    if not isinstance(info, pd.DataFrame):
        print(f'''First argument should be a Pandas DataFrame not a {type(info)}''', file=sys.stderr)
        return None

    curr_dir = os.getcwd()
    os.makedirs(f'{curr_dir}/clusters/hhmdb/aln')
    for x in info.c80e3.unique():
        if info[info.c80e3 == x].c100i100.nunique() >= i:
            os.mkdir(f'{curr_dir}/clusters/{x}')
            os.chdir(f'{curr_dir}/clusters/{x}')
            rdbs.sequence(info[info.c80e3 == x].c100i100.drop_duplicates().to_list(), local_database_path=f'{curr_dir}/{index}').align(method='linsi').to_file(f'{x}.{c80e3}.aln')
        os.chdir(f'{curr_dir}')

    os.chdir('clusters')
    os.system('for x in */*.aln; do y=$(echo $x|cut -f1 -d"/"); echo \#$y|cat - $x > ./hhmdb/aln/$y.aln;done')

    os.chdir('hhmdb')
    os.system(f'/home/leep/ggnicastro/bin/build/build_hhdb.sh ./aln ../{index}')

    os.chdir('aln')
    os.system(f'for x in *.aln; do hhsearch -i $x -d ../{index} -M 50;done')
    #os.system('python3 /home/leep/ggnicastro/bin/hhsearch_table.py')
