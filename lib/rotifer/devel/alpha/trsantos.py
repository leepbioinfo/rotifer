import numpy as np
import pandas as pd
from collections import Counter
from rotifer.db import ncbi
from rotifer.db.ncbi import entrez
from rotifer.devel.alpha import epsoares as rdae
from rotifer.intervals import utils as riu

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

def create_attable(ids, update=False, ndf=None, ipg=None):
    """
    Create an enriched sequence object with attribute table, neighborhoods, and Pfam architecture.

    This function initializes a sequence object, fetches genome neighborhoods, IPG information,
    and taxonomy, and compiles them into a unified attribute table (.df). Neighborhoods are
    mapped only for query proteins, and Pfam domains are integrated from HMM scans.

    Parameters
    ----------
    ids : list
        Accession IDs to build the attribute table from.
    update : bool, optional
        If True, updates the ETE3 taxonomy database before fetching.
    ndf : pd.DataFrame, optional
        Precomputed genome neighborhood DataFrame. If given, skips fetching.
    ipg : pd.DataFrame, optional
        Precomputed IPG information. If given, skips fetching.

    Returns
    -------
    a : rotifer.devel.beta.sequence
        Sequence object with:
        - a.df : main attribute table with assembly, taxonomy ranks, neighborhood info, and Pfam architecture
        - a.ndf : genome neighborhood DataFrame
        - a.ipg : IPG information
        - a.tax : taxonomy table (pivoted)
        - a.hsdf : Pfam HMM scan results
    """
    
    # Step 1: Initialize sequence object
    print("Step 1: Initializing sequence object...")
    a = rdbs.sequence(ids)
    print("Step 1 done.")

    # Step 2: Fetch IPG info if not provided
    if ipg is None:
        print("Step 2: Fetching IPG info...")
        a.ipg = ic.fetchall(ids)
        print("Step 2 done.")
    else:
        a.ipg = ipg
        print("Step 2 skipped (using provided IPG).")

    # Step 3: Fetch genome neighborhood if not provided
    if ndf is None:
        print("Step 3: Fetching genome neighborhood...")
        a.ndf = rdae.add_arch_to_df(gnc.fetchall(ids, ipgs=a.ipg))
        print("Step 3 done.")
    else:
        a.ndf = ndf
        print("Step 3 skipped (using provided NDF).")

    # Step 4: Fetch taxonomy info
    print("Step 4: Fetching taxonomy info...")
    if update:
        print("Updating ETE3 taxonomy database...")
        tc.cursors['ete3'].update_database()
    unique_taxids = a.ndf['taxid'].drop_duplicates()
    tax = tc.fetchall(unique_taxids)
    print("Step 4 done.")

    # Step 5: Expand lineages
    print("Step 5: Expanding lineages...")
    lineage_dict = tc.cursors['ete3'].ete3.get_lineage_translator(unique_taxids.tolist())
    lineage_series = pd.Series(lineage_dict).explode()
    z = lineage_series.reset_index().rename(columns={'index': 'taxid', 0: 'lineage'})
    rank_dict = tc.cursors['ete3'].ete3.get_rank(z['lineage'].tolist())
    taxon_dict = tc.cursors['ete3'].ete3.get_taxid_translator(z['lineage'].tolist())
    z['rank'] = z['lineage'].map(rank_dict)
    z['taxon'] = z['lineage'].map(taxon_dict)
    print("Step 5 done.")

    # Step 6: Pivot taxonomy table
    print("Step 6: Pivoting taxonomy table...")
    ranks = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxon_df = (
        z[z['rank'].isin(ranks)]
        .pivot(index='taxid', columns='rank', values='taxon')
        .reset_index()
    )
    print("Step 6 done.")

    # Step 7: Merge taxonomy to main df
    print("Step 7: Mapping taxonomy and assembly to .df...")
    tmp_tax = a.ndf[['pid', 'taxid', 'assembly']].drop_duplicates()
    tmp_tax = tmp_tax.merge(taxon_df, on='taxid', how='left')
    a.df = a.df.merge(tmp_tax, left_on='id', right_on='pid', how='left', suffixes=('', '_y'))
    a.df.drop(columns=['pid'], inplace=True, errors='ignore')
    a.tax = taxon_df
    print("Step 7 done.")

    # Step 8: Map neighborhoods only for query proteins
    print("Step 8: Mapping neighborhood and ss_neighborhood...")
    cndf = a.ndf.compact()          # default compact
    cndf_ss = a.ndf.compact(strand=True)  # strand-specific

    # Map block_id to .df only for query proteins
    query_block_map = a.ndf[a.ndf['query'] == 1][['pid', 'block_id']].drop_duplicates()
    a.df.loc[a.df['id'].isin(query_block_map['pid']), 'block_id'] = \
        a.df['id'].map(dict(zip(query_block_map['pid'], query_block_map['block_id'])))

    # Map neighborhoods using index-based mapping (clean)
    mask = a.df['block_id'].notna()
    a.df.loc[mask, 'neighborhood'] = a.df.loc[mask, 'block_id'].map(cndf['compact'])
    a.df.loc[mask, 'ss_neighborhood'] = a.df.loc[mask, 'block_id'].map(cndf_ss['compact'])

    # Drop block_id column (optional)
    a.df.drop(columns=['block_id'], inplace=True, errors='ignore')

    # Step 9: Fill missing values
    a.df.fillna('Missing', inplace=True)
    a.tax.fillna('Missing', inplace=True)

    print("Step 8-9 done. Attribute table ready.")

    # Step 10: Pfam integration
    hsdf = rdae.hmmscan(a.df.id.to_list())
    a.hsdf = hsdf
    hsdff = riu.filter_nonoverlapping_regions(hsdf, **riu.config['hmmer'])
    idx = hsdf.groupby('sequence')['evalue'].idxmin()
    best_hits = hsdf.loc[idx, ['sequence', 'model']].set_index('sequence')['model']
    a.df['arch'] = a.df['id'].map(best_hits).fillna('Missing')

    print("Step 10 done. Pfam column added.")

    return a

