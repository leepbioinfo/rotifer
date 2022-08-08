#!/usr/bin/env python3


def dataframe_to_community(df, source='query', target='hit', weight='probability', directed=False):
    '''
    Extract communities from a dataframe using the Louvain algorithm.
    '''
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

    G = nx.from_pandas_edgelist(df[[source, target, weight]], edge_attr=weight)
    partition = community.best_partition(G, weight=weight)
    d = (
        pd.DataFrame.from_dict(partition, orient="index")
        .reset_index()
        .rename({"index": "nodes", 0: "community"}, axis=1)
    )
    return d

