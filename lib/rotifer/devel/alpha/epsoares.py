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

    seqobj.df.sequence = seqobj.df.apply(lambda x: x.sequence[int(x[start])-1:int(x[end])], axis=1)
    seqobj.df['pid'] = seqobj.df['id']
    seqobj.df['id'] = seqobj.df['id'] + '/' + seqobj.df[start].astype(str) + '-' + seqobj.df[end].astype(str)
    seqobj._reset()

    return seqobj

def to_network(df, target=['pfam'], ftype='CDS', interaction=True, ignore = []):
    if isinstance(target, str):
        target = [target]

    w = df.query(f'type == "{ftype}"')[['block_id', 'strand']].copy()
    w['source'] = df[target[0]]

    for col in target[1:]:
        w['source'] = np.where(w['source'].isna(), df[df.type == ftype][col], w['source'])

    w['source'] = w['source'].fillna("?").str.split('+')
    w = w.explode(column='source')
    w['target'] = w['source'].shift(-1)

    if ignore:
        w = w[~(w.source.isin(ignore) | w.target.isin(ignore))].copy()

    w = w[(w.block_id == w.block_id.shift(1)) & (w.strand == w.strand.shift(1))].copy()
    w.loc[w.strand == -1, ['source', 'target']] = w.loc[w.strand == -1, ['target', 'source']].values

    if interaction:
        w['interaction'] = np.where((w.index.to_series() == w.index.to_series().shift(1)), 'fusion', 'neighbor')
        w = w.groupby(['source', 'target', 'interaction'])
    else:
        w = w.groupby(['source', 'target'])

    w = w.agg(weight=('strand', 'count'), blocks=('block_id', 'nunique')).reset_index()
    return w

def compact_for_treeviewer(
        ndf,
        acc,
        columns=['pid','assembly','nucleotide','block_id','organism','lineage','classification','pfam','aravind','compact'],
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
    ndf_acc['compact1'] = ndf_acc.block_id.map(ndfc.compact.to_dict())
    ndf_acc['compacts'] = ndf_acc.block_id.map(ndfcs.compact.to_dict())
    table = ndf_acc[columns].drop_duplicates(columns[0])
    if save:
        table.to_csv(save, sep = '\t', index = False)
        print(f'Table saved to {save}')
    else:
        return table