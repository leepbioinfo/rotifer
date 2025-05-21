import os
import sys
import numpy as np
import pandas as pd

def get_matrix(df, filter_list, rows, columns, n=10, filter_by='pid'):
        import pandas as pd
        filtered_df = df[df[filter_by].isin(filter_list)]
        frequent_rows = filtered_df[rows].value_counts(rows).nlargest(n).index.tolist()
        filtered_df = filtered_df[filtered_df[rows].isin(frequent_rows)]
        matrix = filtered_df.pivot_table(index=rows, columns=columns, aggfunc='size', fill_value=0)
        return matrix

def update_lineage(ndf, preferred="/home/leep/epsoares/projects/databases/data/preferred_taxa.txt"):
    '''
    A function to update the lineage of Gene Neighborhood cursor;
    '''
    from rotifer.db import ncbi
    from rotifer.taxonomy import utils as rtu
    tc = ncbi.TaxonomyCursor()
    tax = tc.fetchall(ndf.taxid.dropna().drop_duplicates().tolist())
    tax['taxid'] = tax.taxid.astype(int)
    ndf.classification = ndf.taxid.map(tax.set_index('taxid').classification.to_dict())
    ndf.lineage = rtu.lineage(ndf.classification, preferred_taxa=[ x.strip() for x in open(preferred,"rt")]).tolist()
    return ndf

def extract_envelope(df, seqobj=None, start='estart', end='eend', expand=10, local_database_path='/databases/fadb/nr/nr'):
    '''
    Function to extract the envelope of a model match.
    '''

    import pandas as pd
    import rotifer.devel.beta.sequence as rdbs
    import numpy as np

    if len(df) == 0:
        print(f'Empty input')
        seqobj = rdbs.sequence()

    if seqobj == None:
        seqobj = rdbs.sequence(df.sequence.drop_duplicates().to_list(), local_database_path=local_database_path)

    seqobj.df = seqobj.df.merge(df.rename({'sequence':'id'},axis=1), how='left')
    seqobj.df.end = seqobj.df[end].fillna(seqobj.df.length)
    seqobj.df.start = seqobj.df[start].fillna(1)

    if expand:
        seqobj.df['start'] = np.where(seqobj.df.start<=expand,1,seqobj.df.start-expand+1)
        seqobj.df['end'] = np.where(seqobj.df.end+expand > seqobj.df.length, seqobj.df.length, seqobj.df.end+expand)

    seqobj.df.sequence = seqobj.df.apply(lambda row: row['sequence'][row['start'] - 1:row['end']], axis=1)
    seqobj.df['pid'] = seqobj.df['id']
    seqobj.df['id'] = seqobj.df['id'] + '/' + seqobj.df[start].astype(str) + '-' + seqobj.df[end].astype(str)
    seqobj._reset()

    return seqobj

def to_network(df, target=['pfam'], ftype=['CDS'], interaction=True, ignore = [], strand = True):
    if isinstance(ftype, str):
        ftype = [ftype]
    if isinstance(target, str):
        target = [target]

    w = df.copy()
    if strand:
        w = df.neighbors(strand='same')
        w['rid'] = list(range(1,len(w)+1))
        w.rid = w.rid * w.strand
        w.sort_values(['rid'], inplace=True)

    # Building the source column
    w = w.query(f'type == "{ftype}"').block_id.reset_index().drop('index', axis=1)
    w['source'] = df[target[0]]
    for col in target[1:]:
        w['source'] = np.where(w['source'].isna(), df[df.type == ftype][col], w['source'])
    w['source'] = w['source'].fillna("?").str.split('+')
    w = w.explode(column='source')
    if ignore:
        w = w[~w.source.isin(ignore)].copy()

    # Building target data
    w['tblock_id'] = w['block_id'].shift(-1)
    w['target'] = w['source'].shift(-1)

    # Selecting same block rows
    w = w[w.block_id == w.tblock_id.shift(1)].copy()

    # Fix source target order when not restricted to the same strand
    sameprotein = (w.index.to_series() == w.index.to_series().shift(1))
    if not strand:
        reverse = (~sameprotein) & (w.source > w.target)
        w.loc[reverse,['source','target']] = w.loc[reverse, ['target', 'source']].values

    if interaction:
        w['interaction'] = np.where(sameprotein, 'fusion', 'neighbor')
        w = w.groupby(['source', 'target', 'interaction'])
    else:
        w = w.groupby(['source', 'target'])

    w = w.agg(weight=('block_id', 'count'), blocks=('block_id', 'nunique')).reset_index()
    return w

def compact_for_treeviewer(
        ndf,
        acc,
        columns=['pid','assembly','nucleotide','block_id','organism','lineage','classification','pfam','aravind','compact','compacts'],
        save=None,
    ):
    '''
    Add a compact GeneNeighborhoodDF representation, select
    and reorder columns to match those required by TreeViewer
    and FigTree (leaf identifier as first column).
    '''
    ndfc = ndf.compact()
    ndfcs = ndf.select_neighbors(strand = True).compact()
    ndf_acc = ndf[ndf.pid.isin(acc)]
    ndf_acc['compact'] = ndf_acc.block_id.map(ndfc['compact'].to_dict())
    ndf_acc['compacts'] = ndf_acc.block_id.map(ndfcs['compact'].to_dict())
    table = ndf_acc[columns].drop_duplicates(columns[0])
    if save:
        table.to_csv(save, sep = '\t', index = False)
        print(f'Table saved to {save}')
    else:
        return table

def make_heatmap(
        ndf,
        name = None,
        domain_list = [],
        format_table = True,
        color_list = ['white','blue','green','red'],
        cbar = False,
        fmt = 'd',
        annot = False,
        linewidths = 1,
        tree_file = None,
        tree = False,
        figsize=None
    ):
    '''
    Creates a heatmap from a dataframe output of full annotate
    to represent main occurrence. The colors are correspondent 
    to the values in the dataframe, by exemple, four colors 
    corresponds to the values: 0, 1, 2 and 3, respectively. 
    The columns in the dataframe will be the columns in the 
    heatmap, and each line represents a occurrence from the 
    respective domain and the colors the type of correspondence, 
    that is determined by the numbers utilized. 0 equals to non 
    correlation, 1 correlation by neighborhood, 2 correlation by 
    fusions and 3 by both. If a tree file is delivered by the user
    the heatmap will sort the occurrences by the order of the tree.
    '''
 
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import pandas as pd
              
    if format_table == True:
        hm = ndf.set_index('pid')[['block_id']]
        hm['neighbors'] = hm.block_id.map(ndf.query('query == 0').groupby('block_id').agg(neighbors = ('pfam','sum')).neighbors.to_dict())
        hm['pfam'] = hm.block_id.map(ndf.query('query == 1').set_index('block_id').pfam.to_dict())
        for x in domain_list:
            inpfam = hm['pfam'].str.contains(x, na = False)
            inneighbors = hm['neighbors'].str.contains(x, na = False)
            hm[f'{x}'] = ((inpfam & inneighbors).astype(int)*3 + (inpfam & ~inneighbors).astype(int)*2 + (~inpfam & inneighbors).astype(int)*1 + (~inpfam & ~inneighbors).astype(int)*0)
        df = hm[domain_list]
                                                                   
    if tree == True:
        import ete3
        t = ete3.Tree(f'{tree_file}')
        df = df.reindex(df.reindex(t.get_leaf_names()).index.str.replace('\'','')).fillna(0).astype(int)
    
    if figsize:
        figsize=figsize
    else:
        figsize=(len(df.columns),len(df)/5)
  
    cmap = colors.ListedColormap(color_list)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(data = df, cmap = cmap, cbar = cbar, fmt = fmt, annot = annot, linewidths = linewidths, ax=ax)
    plt.savefig(f'{name}', bbox_inches='tight')
    plt.close(fig)

