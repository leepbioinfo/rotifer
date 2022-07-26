#!/usr/bin/env python3

import os
import sys
import rotifer
import pandas as pd
from rotifer.db.ncbi import NcbiConfig
logger = rotifer.logging.getLogger(__name__)

# Load Identical Protein Reports
def ipg(ncbi, batch_size=200, *args, **kwargs):
    '''
    Retrieve and annotate NCBI's Identical Protein Reports (IPG).

    Usage:
      a = ncbi()
      a.submit(['WP_028373859.1','WP_043290818.1','WP_082188481.1'])
      b = a.read(method='ipg')

    Returns:
      Pandas DataFrame

    Parameters:
      batch_size : integer
        Number of IPGs to retrieve at each batch
    '''

    # Method dependencies
    import numpy as np
    import pandas as pd
    from Bio import Entrez
    Entrez.email = ncbi.email
    cols = ['id','ipg_source','nucleotide','start','stop','strand','pid','description','ipg_organism','strain','assembly']
    added = ['order','is_query','representative']
    emptyDF = pd.DataFrame(columns=cols+added)
    logger.info(f'downloading IPG reports for {len(ncbi)} protein accessions...')

    # Backup and check queries
    queries = ncbi.submit()
    if len(queries) == 0:
        return emptyDF
    logger.info(f'searching {len(ncbi)} accessions...')

    # Using Efetch directly with batch_size IDs per batch!
    n = 1
    ipgs = []
    pos = list(range(0,len(queries),batch_size))
    for s in pos:
        e = s + batch_size
        if e > len(queries):
            e = len(queries)
        batch = queries[s:e]
        logger.info(f'downloading batch {n} ([{s}:{e}]) of {len(pos)}')
        handle = None
        try:
            handle = Entrez.efetch(db='protein', rettype='ipg', retmode='text', api_key=NcbiConfig['api_key'], id = ",".join(batch))
        except RuntimeError:
            logger.info(f'batch {n}, runtime error: '+str(sys.exc_info()[1]))
            continue
        except:
            logger.info(f'batch {n}, exception: '+str(sys.exc_info()[1]))
            continue
        if handle:
            try:
                ipg = pd.read_csv(handle, sep='\t', names=cols, header=0).drop_duplicates()
            except:
                logger.info(f'batch {n}, could not read IPG: '+str(sys.exc_info()[1]))
                continue
            handle.close()
            if not isinstance(ipg,pd.DataFrame) or ipg.empty:
                logger.info(f'batch {n} has no IPGs! Download or parsing error? Moving to next batch...')
            numeric = pd.to_numeric(ipg.id, errors="coerce")
            errors = numeric.isna()
            if errors.any():
                logger.info(f'Errors in IPG for batch {n}:\n'+ipg[errors].to_string())
                continue
            ipg.id = numeric
            ipg = ipg[~errors]
            if ipg.empty:
                logger.info(f'after removing errors, batch {n} was found to have no IPGs! Ignoring...')
                continue
            o = pd.Series(range(1, len(ipg) + 1))
            c = pd.Series(np.where(ipg.id != ipg.id.shift(1), o.values, pd.NA)).ffill()
            ipg['order'] = (o - c).values
            ipgs.append(ipg)
        n += 1

    # Failure! Give up: no IPG was downloaded! Database failure?
    if not ipgs:
        logger.info(f'{len(ncbi)} queries but 0 IPGs! Make sure your queries belong to the NCBI protein database.')
        return emptyDF

    # Success! Concatenate all batches
    ipgs = pd.concat(ipgs).reset_index(drop=True).sort_values(['id','order']).drop_duplicates()
    found = set(queries).intersection(set(ipgs.pid.unique()))
    missing = set(queries) - set(found)
    logger.info(f'{len(ipgs)} IPG rows fetched, {len(found)} queries found, {len(missing)} missing accessions.')

    # Filter error messages
    if ipgs.empty:
        ncbi.missing(missing)
        return emptyDF

    # Register query proteins
    ipgs['is_query'] = ipgs.pid.isin(queries).astype(int)

    # Process batches of different sizes
    # More than one query?
    if len(queries) > 1: # Many queries
        #  Register first query protein as representative
        rep = ipgs.loc[ipgs['is_query'] == 1,['id','pid']].drop_duplicates('id', keep='first')
        idmap = {k:v for k,v in zip(rep['id'].values, rep['pid'].values) }
        ipgs['representative'] = ipgs['id'].map(idmap)

        # Remove all IPGs that have no representative
        ipgs = ipgs[ipgs['representative'].notnull()]

        # Multiple queries and still missing some? Let's try to fetch one at a time...
        if len(missing) > 0:
            for lost in missing:
                ncbi.submit(lost)
                ipg = ncbi.read('ipg')
                if ipgs.empty:
                    ipgs = ipg
                else:
                    ipgs = ipgs.append(ipg, ignore_index=True).drop_duplicates().reset_index(drop=True)
            ncbi.submit(queries) # Reset query list

    # One query
    else:
        ipgs['representative'] = queries[0]

    # Set id to numeric and update list of missing queries
    ipgs['id'] = pd.to_numeric(ipgs.id)
    missing = set([x for x in queries if x not in ipgs.pid.append(ipgs.representative).unique()])
    ncbi.missing(missing)
    found = set([x for x in queries if x not in missing])
    logger.info(f'{len(ipgs)} IPG rows fetched, {len(found)} queries found, {len(missing)} missing accessions.')

    # Return the expected dataframe
    return ipgs

# Load taxonomy data as a simple dataframe (tree as linearized path)
def taxonomy(ncbi, fetch=['ete3','entrez'], missing=False, ete3=None, preferred_taxa=None, *args, **kwargs):
    '''
    Load NCBI taxonomy data

    Usage:
      a = ncbi()
      a.submit('9606')
      b = a.read('taxonomy')

      or, by organism name using pre-loaded ete3 object,

      from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa
      e = NCBITaxa()
      a = ncbi()
      a.submit(['Homo sapiens'])
      b = a.read('taxonomy', fetch=['ete3'], ete3=e, query_type='scientific')

    Returns:
      Pandas DataFrame
      Columns: taxid, organism, superkingdom, lineage, taxonomy

    Parameters:
      fetch    : list of fetchers
      missing  : retrieve missing entries only
      ete3     : ete3's NCBITaxa object
                 If set to true, a new NCBITaxa object is created
    '''

    # Set log format
    logger.info(f'main: loading assembly reports...')

    # Loop over methods
    stack = []
    queries = ncbi.submit().copy()
    if missing:
        missing = ncbi.missing()
    else:
        missing = queries.copy()
    for method in fetch:
        if not missing:
            break
        ncbi.submit(missing)
        if method == 'ete3':
            tax = __taxonomy_from_ete3(ncbi, ete3=ete3, preferred_taxa=preferred_taxa)
        elif method == 'entrez':
            continue
        missing = ncbi.missing().copy()
        stack.append(tax)
        logger.info(f'{len(tax)} taxa loaded by {method}.')
    if len(stack) > 0:
        tax = pd.concat(stack).drop_duplicates()
    else:
        tax = pd.DataFrame(columns='taxid organism superkingdom lineage taxonomy'.split(' '))

    # Return
    ncbi.submit(queries)
    ncbi.missing(missing)
    logger.info(f'main: {len(tax.taxid.unique())} taxids found and {len(missing)} missing!')
    return tax

# Internal methods

def __taxonomy_from_ete3(ncbi, ete3=None, preferred_taxa=None):
    # Make sure ete3 is a ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa
    if not isinstance(ete3,NCBITaxa):
        ete3 = NCBITaxa()
    logger.info(f'Searching for {len(ncbi.submit())} taxa.')

    # Fetch lineages
    i = 0
    l = set() # Non-redundant list of taxa
    s = {}    # Dictionary: query to lineage map
    missing = []
    query = ncbi.submit().copy()
    nq = len(query)
    while i < nq:
        try:
            query[i] = int(query[i])
        except:
            pass
        if not isinstance(query[i],int):
            d = ete3.get_name_translator([query[i]])
            if query[i] not in d:
                missing.append(query[i])
                i = i + 1
                continue
            query.extend([ x for x in d[query[i]][1:] if x not in query ])
            query[i] = d[query[i]][0]
            nq = len(query)
        try:
            v = ete3.get_lineage(query[i])
        except:
            missing.append(query[i])
            i = i + 1
            continue
        l.update(v)
        s[query[i]] = v
        i = i + 1
    logger.info(f'loaded {len(s)} lineages.')

    # Fetch names
    l = ete3.get_taxid_translator(list(l))
    logger.info(f'loaded {len(l)} taxon names.')

    # Translate all lineages
    cols = ['taxid','organism','superkingdom','taxonomy']
    data = { k: [] for k in cols }
    for x in query:
        # Register missing taxids
        if x not in s:
            #missing.append(x)
            continue

        # Translate each taxa
        names = [ l[v] for v in s[x] ]
        if names[0] == 'root':
            del(names[0])
        if names[0] == 'cellular organisms':
            del(names[0])
        full = "; ".join(names)
        r = [x, names[-1], names[0], full]

        # Store data
        for j in list(range(4)):
            data[cols[j]].append(r[j])

    # Let the user knows we got here
    logger.info(f'translated {len(data[cols[0]])} lineages.')

    # Cleanup: remove missing taxids and register
    data = pd.DataFrame(data)
    from rotifer.taxonomy.utils import lineage
    data.insert(3, "lineage", lineage(data.taxonomy, preferred_taxa=preferred_taxa))
    ncbi.missing(missing)
    return data

# Is this library being used as a script?
if __name__ == '__main__':
    pass
