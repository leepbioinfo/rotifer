import numpy as np
import pandas as pd
from collections import Counter
from rotifer.db import ncbi
from rotifer.db.ncbi import entrez

def taxon_summary(
    df,
    rank=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
    update=False,
    new_ndf=False
):
    """
    Generate a taxonomic summary from a DataFrame with NCBI TaxIDs, including query and assembly counts
    grouped by the specified taxonomic rank.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame (usually a ndf) containing at least a 'taxid' column (NCBI Taxonomy ID),
        a 'query' column (1 or 0 indicating query status), and an 'assembly' column (assembly accession).
        
    rank : list of str, optional
        List of NCBI taxonomic ranks to extract. The final rank in the list (e.g. 'species') will be used
        to group query and assembly counts. Defaults to the standard 7-rank Linnaean hierarchy:
        ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'].
        
    update : bool, optional
        If True, updates the local NCBI taxonomy database via ETE3 before fetching lineages.
        This is recommended the first time you run the function or if your taxonomy data may be outdated.

    new_ndf : bool, optional
        If True, returns the original DataFrame with taxonomic columns (from 'rank') merged in.
        Not recommended to run this function on the new_ndf created, usually results in error.

    Returns
    -------
    a : pandas.DataFrame
        The original DataFrame with taxonomic columns (from `rank`) merged in.

    taxon_df : pandas.DataFrame
        A summary DataFrame, one row per unique taxon at the final rank level,
        with columns:
        - 'taxid' and taxonomic ranks (as defined in `rank`)
        - 'query_tax_count': Number of query sequences mapped to each taxon
        - 'query_tax_pct': Percentage of queries mapped to each taxon
        - 'tax_assembly_count': Number of unique assemblies mapped to each taxon
        - 'tax_assembly_pct': Percentage of assemblies mapped to each taxon
    """
    
    tc = ncbi.TaxonomyCursor() 
    
    # Update database so we can use the most up-to-date taxonomy
    if update:
        tc.cursors['ete3'].update_database()
    
    a = df.copy() 

    # Fetching all taxids
    ndftax = tc.fetchall(a.taxid.drop_duplicates().tolist()) 
    
    # Creating a table to display the requested taxonomies
    z = pd.Series(tc.cursors['ete3'].ete3.get_lineage_translator(ndftax.taxid.drop_duplicates().tolist()), name='lineage').reset_index().rename({'index':'taxid'}, axis=1).explode('lineage')
    z['rank'] = z.lineage.map(tc.cursors['ete3'].ete3.get_rank(z.lineage.drop_duplicates().tolist()))
    z['taxon'] = z.lineage.map(tc.cursors['ete3'].ete3.get_taxid_translator(z.lineage.drop_duplicates().tolist()))
    taxon_df = z[z['rank'].isin(rank)].pivot(index='taxid', columns='rank', values='taxon').reset_index()
    
    # Merging the summary to ndf
    a = a.merge(taxon_df, on='taxid', how='left')
    
    # Measuring the amount of queries that belong to each taxa
    query_df = a[a['query'] == 1]
    query_counts = (query_df.groupby(rank[-1]).size().rename('query_tax_count').reset_index())
    total_queries = query_counts['query_tax_count'].sum()
    query_counts['query_tax_pct'] = (query_counts['query_tax_count'] / total_queries) * 100

    # Measuring the amount of assemblies that belong to each taxa
    assembly_df = a[['assembly', rank[-1]]].drop_duplicates()
    assembly_counts = (assembly_df.groupby(rank[-1]).size().rename('tax_assembly_count').reset_index())
    total_assemblies = assembly_counts['tax_assembly_count'].sum()
    assembly_counts['tax_assembly_pct'] = (assembly_counts['tax_assembly_count'] / total_assemblies) * 100

    taxon_df = taxon_df.merge(query_counts, on=rank[-1], how='left')
    taxon_df = taxon_df.merge(assembly_counts, on=rank[-1], how='left')

    taxon_df[['query_tax_count', 'query_tax_pct', 'tax_assembly_count', 'tax_assembly_pct']] = (taxon_df[['query_tax_count', 'query_tax_pct', 'tax_assembly_count', 'tax_assembly_pct']].fillna(0))
    taxon_df['query_tax_count'] = taxon_df['query_tax_count'].astype(int)
    taxon_df['tax_assembly_count'] = taxon_df['tax_assembly_count'].astype(int)    
    
    taxon_df = taxon_df.drop_duplicates(subset=[rank[-1]])
    columns = ['taxid'] + rank + ['query_tax_count', 'query_tax_pct', 'tax_assembly_count', 'tax_assembly_pct']
    taxon_df = taxon_df[columns]
    taxon_df = taxon_df.sort_values(by=['tax_assembly_pct', 'query_tax_pct'], ascending=[False, False])
    
    if new_ndf:
        return a, taxon_df
    else:
        return taxon_df

def shannon(self, ignore_gaps=True):
    """
    Computes Shannon entropy for each column in the alignment.

    Parameters
    ----------
    ignore_gaps : bool
        Whether to exclude gaps ('-' or '.') from entropy calculation.

    Notes
    -----
    You need to run the following commands in order for this method to work properly:

    from rotifer.devel.beta import sequence as rdbs
    rdbs.sequence.compute_shannon_entropy = compute_shannon_entropy

    Returns
    -------
    entropy : pd.Series
        A pandas Series with entropy values indexed by column number.
    """

    # Make sure we only work with sequences
    sequences = self.df[self.df['type'] == 'sequence']['sequence'].tolist()

    if not sequences:
        raise ValueError("No sequences found to compute entropy.")

    # Transpose the alignment into columns
    alignment_length = len(sequences[0])
    columns = zip(*sequences)  # Converts list of sequences to columns

    entropy_values = []
    for col in columns:
        if ignore_gaps:
            col = [res for res in col if res not in ['-', '.']]
        freqs = Counter(col)
        total = sum(freqs.values())
        if total == 0:
            entropy = 0
        else:
            probs = [count / total for count in freqs.values()]
            entropy = -sum(p * np.log2(p) for p in probs if p > 0)
        entropy_values.append(entropy)

    return pd.Series(entropy_values, name='shannon_entropy')

def flag_best_id(group):
    """
    Flag one "best" id per group.

    Priority: pick a random row with id_type 'RefSeq', else 'EMBL-CDS'.
    Returns the same group DataFrame with a new integer 'flag' column (0/1).
    Example:
        radsamorg = radsamorg.groupby('qid', group_keys=False).apply(flag_best_id)
    """
    n = len(group)
    flags = np.zeros(n, dtype=int)

    id_types = group['id_type'].to_numpy()
    for t in ("RefSeq", "EMBL-CDS"):
        pos = np.flatnonzero(id_types == t)   # positions within the group (0..n-1)
        if pos.size:
            choice = np.random.choice(pos)
            flags[choice] = 1
            break

    group = group.copy()     # avoid SettingWithCopyWarning
    group['flag'] = flags
    return group
