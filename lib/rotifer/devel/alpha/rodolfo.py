#!/usr/bin/env python3


#todo

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

