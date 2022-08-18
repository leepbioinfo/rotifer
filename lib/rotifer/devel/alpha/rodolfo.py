#!/usr/bin/env python3
import rotifer
from rotifer import GlobalConfig
logger = rotifer.logging.getLogger(__name__)

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

def build_hhuite_database(
    df,
    sequences,
    group="c80e3",
    members="c80i70",
    output_directory="hhmdb",
    alignment_method='famsa',
    minimum_group_size=3,
    threads=10,
):
    '''
    Build a HH-suite database for all protein groups.

    Parameters:
      df: a Pandas Dataframe
        The members column should list all sequences the user
        wishes to add to the alignment of each group
      sequences: string
        Path to a FASTA file with the sequences under analysis.
      group: string, default 'c80e3'
        Name of the column containing group identifiers
      members: string, default 'c80i70'
        Name of the column storing sequence identifiers
      output_directory: path-list string, default 'hhmdb'
        Name of the directory to save all output files
      minimum_group_size: integer, default 3
        Size of the smallest cluster that must be added
        to the database
      threads: integer, default 10
        Number os CPUs to use while building alignments
    '''
    import os
    import tempfile
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence

    # Create output directory
    if not os.path.exists(f'{output_directory}/aln'):
        os.makedirs(f'{output_directory}/aln')

    # Build MSAs
    for group_cluster in df[group].unique():
        g = df.loc[df[group] == group_cluster,members].drop_duplicates().dropna()
        if g.nunique() >= minimum_group_size:
            aln = sequence(g.tolist(),"accession", name=group_cluster, local_database_path=sequences)
            aln = aln.align(method="famsa")
            aln.to_file(f'{output_directory}/aln/{group_cluster}.a3m','a3m')

    # Build hh-suite database
    with tempfile.TemporaryDirectory() as td:
        cmd = rotifer.GlobalConfig['base'] + f'/share/rotifer/scripts/build_hhdb.sh'
        Popen(f'{cmd} {td}/{output_directory} {output_directory}/aln').communicate()
        for orig in ['cs219','a3m.ordered','hhm.ordered']:
            dest = orig.replace(".ordered","")
            os.rename(f'{td}/{output_directory}.{orig}.ffdata',f'./{output_directory}/{output_directory}_{dest}.ffdata')
            os.rename(f'{td}/{output_directory}.{orig}.ffindex',f'./{output_directory}/{output_directory}_{dest}.ffindex')

def hhblits_all(input_directory, output_directory=None, input_suffix="a3m", output_suffix="hhr"):
    """
    Run all alignments in a directory against a HH-suite database.
    """
    from glob import glob
    if not output_directory:
        output_directory = f'{input_directory}/{output_suffix}'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    for aln in glob(f'{input_directory}/*.{input_suffix}'):
        y = os.path.basename(aln).replace(".{input_suffix}","")
        cmd = f'hhsearch --cpu {cpu} -i {aln} -d {input_directory}/{input_directory} -o {output_directory}/{y}.{output_suffix}'
        Popen(cmd, stdout=PIPE, shell=True).communicate()

