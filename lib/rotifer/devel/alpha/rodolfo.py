#!/usr/bin/env python3

def annotate_community(df, column='community', reference='query', ):
    '''
    add or modify communities in a dataframe
    Parameters:
        df: Pandas dataframe
        dataframe to be annotated
        column:
        rows: which rows to annotate
        reference: 
    '''
    # loc to change info
    # or map
    # df.loc[rows, column] = df[reference].map(dict)
    # dict is a result of a df2comm
    # cut communities by cluster size and inputs two tables
    # annotated table and graph
    return


def dataframe_to_community(df, source='query', target='hit', weight='probability', directed=False):
    '''
    Extract communities from a dataframe using the Louvain algorithm.
    '''
    import pandas as pd
    import numpy as np
    import networkx as nx
    import community

    # The following lines are responsible for transforming the input matrix
    # into a symmetric matrix, i.e. an undirected graph
    if not directed:
        df[source], df[target] = np.where(
            df[source] > df[target],
            [df[source], df[target]],
            [df[target], df[source]],
        )
        df = df.sort_values(weight).drop_duplicates([source, target])

    G = nx.from_pandas_edgelist(df.filter([source, target, weight]), source=source, target=target, edge_attr=weight)
    partition = community.best_partition(G, weight=weight)
    d = (
        pd.DataFrame.from_dict(partition, orient="index")
        .reset_index()
        .rename({"index": "nodes", 0: "community"}, axis=1)
    )
    return d

def hhmdb(
    df,
    esl_index_file,
    grouper="c80e3",
    redundancy_cluster="c80i70",
    cpu=10,
    prefix="hhmdb",
    size=3,
):
    '''
    Builds a3m and hhm database for the chosen cluster in the current directory
    '''

    import tempfile
    from subprocess import Popen, PIPE, STDOUT
    import os
    from rotifer.devel.beta.sequence import sequence

    curr_dir = os.getcwd()
    if not os.path.exists(curr_dir+"/hhmdb/aln/"):
        os.makedirs(curr_dir+"/hhmdb/aln/")

    with tempfile.TemporaryDirectory() as tdf:
        for group_cluster in df[grouper].unique():
            if df[df[grouper] == group_cluster][redundancy_cluster].drop_duplicates().dropna().nunique() >= size :
                df[df[grouper] == group_cluster][
                redundancy_cluster
                ].drop_duplicates().dropna().to_csv(
                f"{tdf}/accs", index=None, header=None
                )
                Popen(
                f"esl-sfetch -f {esl_index_file} {tdf}/accs | mafft --anysymbol --maxiterate 1000 --localpair --thread {cpu} - > {tdf}/accs.fa",
                stdout=PIPE,
                shell=True,
                ).communicate()
                sequence(f"{tdf}/accs.fa").to_file(f"./hhmdb/aln/{group_cluster}.fa")
    Popen(f'for x in ./hhmdb/aln/*.fa; do y=$(echo $x|cut -f 4 -d "/" | sed "s/\.fa//"); echo \#$y|cat - $x > ./hhmdb/aln/$y.aln; rm $x;done', stdout=PIPE, shell=True).communicate()
    with tempfile.TemporaryDirectory() as td:
        Popen(f"ffindex_build -s {td}/{prefix}.msa.ffdata {td}/{prefix}.msa.ffindex ./hhmdb/aln/", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_apply {td}/{prefix}.msa.ffdata {td}/{prefix}.msa.ffindex -i {td}/{prefix}.a3m.ffindex -d {td}/{prefix}.a3m.ffdata -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_apply {td}/{prefix}.a3m.ffdata {td}/{prefix}.a3m.ffindex -i {td}/{prefix}.hhm.ffindex -d {td}/{prefix}.hhm.ffdata -- hhmake -i stdin -o stdout -v 0", stdout=PIPE, shell=True).communicate()
        Popen(f"cstranslate -f -x 0.3 -c 4 -I a3m -i {td}/{prefix}.a3m -o {td}/{prefix}.cs219", stdout=PIPE, shell=True).communicate()
        Popen(f" sort -k3 -n -r {td}/{prefix}.cs219.ffindex | cut -f1 > {td}/{prefix}.sorting.dat", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_order {td}/{prefix}.sorting.dat {td}/{prefix}.hhm.ffdata {td}/{prefix}.hhm.ffindex {td}/{prefix}.hhm.ordered.ffdata {td}/{prefix}.hhm.ordered.ffindex ", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_order {td}/{prefix}.sorting.dat {td}/{prefix}.a3m.ffdata {td}/{prefix}.a3m.ffindex {td}/{prefix}.a3m.ordered.ffdata {td}/{prefix}.a3m.ordered.ffindex ", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.cs219.ffdata ./hhmdb/{prefix}_cs219.ffdata", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.a3m.ordered.ffdata ./hhmdb/{prefix}_a3m.ffdata", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.hhm.ordered.ffdata ./hhmdb/{prefix}_hhm.ffdata", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.cs219.ffindex ./hhmdb/{prefix}_cs219.ffindex", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.a3m.ordered.ffindex ./hhmdb/{prefix}_a3m.ffindex", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.hhm.ordered.ffindex ./hhmdb/{prefix}_hhm.ffindex", stdout=PIPE, shell=True).communicate()
    Popen(f'for x in ./hhmdb/aln/*.aln; do y=$(echo $x|cut -f 4 -d "/" | sed "s/\.aln//"); hhsearch --cpu {cpu} -i $x -d ./hhmdb/{prefix} -o ./hhmdb/aln/$y.hhr;done', stdout=PIPE, shell=True).communicate()

