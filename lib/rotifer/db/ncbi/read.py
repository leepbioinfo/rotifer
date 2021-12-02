#!/usr/bin/env python3

import os
import sys
import logging
import pandas as pd
from rotifer.db.ncbi import NcbiConfig

# Load NCBI assembly reports
def assembly_reports(ncbi, baseurl=f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS', columns=[], query_type='assembly', taxonomy=None, verbose=0, *args, **kwargs):
    '''
    Load NCBI assembly reports from a local directory or FTP.

    Usage:
      a = ncbi()
      b = a.read(method='assembly_reports')

      or, with filters,

      a = ncbi()
      a.submit(['GCF_900504695.1', 'GCF_004636045.1'])
      b = a.read(method='assembly_reports', columns=['id','pid'])

    Returns:
      Pandas DataFrame

    Parameters:
      baseurl    : URL or directory with assembly_summary_*.txt files
      columns    : (optional) list of columns to retrieve
      query_type : column name to use for filtering returned rows
                   Set to None to use pd.DataFrame.query
      taxonomy   : ete3's NCBITaxa object
                   If set to true, a new NCBITaxa object is created

    Extra attributes:
      loaded_from : data source (same as baseurl)
    '''

    # Method dependencies
    import pandas as pd
    from glob import glob

    # Set log format
    logger = logging.getLogger('rotifer.db.ncbi')
    if verbose:
        logger.setLevel(verbose)
        logger.info(f'main: loading assembly reports...')

    # Load assembly reports
    assemblies = list()
    for x in ['refseq', 'genbank', 'refseq_historical', 'genbank_historical']:
        if os.path.exists(baseurl): # Local file
            url = os.path.join(baseurl, f'assembly_summary_{x}.txt')
            if not os.path.exists(url):
                if verbose:
                    logger.warning(f'{__name__}: {url} not found. Ignoring...')
                continue
        else: # FTP
            url = f'{baseurl}/assembly_summary_{x}.txt'
        _ = pd.read_csv(url, sep ="\t", skiprows=[0])
        _.rename({'# assembly_accession':'assembly'}, axis=1, inplace=True)
        _['source'] = x
        _['loaded_from'] = url
        assemblies.append(_)
        if verbose:
            logger.info(f'{__name__}: {url}, {len(_)} rows, {len(assemblies)} loaded')
    assemblies = pd.concat(assemblies, ignore_index=True)
    if verbose:
        print(f'{__name__}: loaded {len(assemblies)} assembly summaries.', file=sys.stderr)

    # Make sure the ftp_path columns refers to the ftp site as we expect
    if 'ftp_path' in assemblies.columns:
        assemblies.ftp_path = assemblies.ftp_path.str.replace('https','ftp')

    # Filter rows and columns in the assemblies dataframe
    queries = ncbi.submit()
    if queries:
        if query_type in assemblies.columns: # Exact match values in a column
            assemblies = assemblies[assemblies[query_type].isin(queries)]
            if verbose:
                print(f'{__name__}: {len(queries)} exact matches required for column {query_type}: {len(assemblies)} rows left.', file=sys.stderr)
        else: # Queries are logical statements
            for query in queries:
                assemblies = assemblies.query(query)
            if verbose:
                print(f'{__name__}: {len(queries)} filters applied: {len(assemblies)} rows left.', file=sys.stderr)

    # Add taxonomy
    if isinstance(taxonomy,pd.DataFrame) or taxonomy:
        if not isinstance(taxonomy,pd.DataFrame):
            ncbi.submit(list(assemblies.taxid.unique()))
            taxonomy = ncbi.read('taxonomy', ete3=taxonomy, verbose=verbose)
        if isinstance(taxonomy,pd.DataFrame):
            assemblies = assemblies.merge(taxonomy, left_on='taxid', right_on='taxid', how='left')
        if verbose:
            print(f'{__name__}: {len(assemblies)} assemblies left-merged with taxonomy dataframe.', file=sys.stderr)

    # Filter columns and return pandas object
    if columns:
        assemblies = assemblies.filter(columns)

    # Reset ncbi object, update missing list and return
    ncbi.submit(queries)
    found = set(list(assemblies.assembly))
    ncbi.missing([ x for x in queries if x not in found ])
    if verbose:
        logger.info(f'main: {len(assemblies)} assembly reports loaded!')
    return assemblies

# Load Identical Protein Reports
def ipg(ncbi, fetch=['entrez'], assembly_reports=False, verbose=False, batch_size=200, *args, **kwargs):
    '''
    Retrieve and annotate NCBI's Identical Protein Reports (IPG).

    Usage:
      a = ncbi()
      a.submit(['WP_028373859.1','WP_043290818.1','WP_082188481.1'])
      b = a.read(method='ipg')

    Returns:
      Pandas DataFrame

    Parameters:
      fetch            : how to download data (see rotifer.db.ncbi.fetch)
      assembly_reports : assembly reports dataframe to set genome priority
                         If set to True, the data is downloaded
    '''

    # Method dependencies
    import pandas as pd
    from Bio import Entrez
    Entrez.email = ncbi.email
    cols = ['id','source','nucleotide','start','stop','strand','pid','description','organism','strain','assembly']
    ipgs = pd.DataFrame(columns=cols)

    # Set log format
    logger = logging.getLogger(__name__)
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.info(f'main: downloading IPG reports for {len(ncbi)} protein accessions...')

    # Backup and check queries
    queries = ncbi.submit()
    if len(queries) == 0:
        return pd.DataFrame(columns=[*cols, 'is_query', 'representative'])
    if verbose:
        print(f'{__name__}: searching {len(ncbi)} accessions...', file=sys.stderr)

    # Using Efetch directly with 200 IDs per batch!
    ipgs = []
    pos = list(range(0,len(queries),batch_size))
    for s in pos:
        # Prepare and submit current batch
        e = s + batch_size
        batch = queries[s:e]
        handle = None
        try:
            handle = Entrez.efetch(db='protein', rettype='ipg', retmode='text', api_key=ncbi.api_key(), id = ",".join(batch))
        except RuntimeError:
            if verbose > 1:
                print(f'Efetch: {len(queries)} queries, {len(missing)} accessions. Runtime error: '+str(sys.exc_info()[1]), file=sys.stderr)
        except:
            if verbose > 1:
                print(f'Efetch: {len(queries)} queries, {len(missing)} accessions. Unexpected error: '+str(sys.exc_info()[0:2]), file=sys.stderr)
        # Merge Efetch results with Epost+Efetch results
        if handle:
            ipg = pd.read_csv(handle, sep='\t', names=cols, header=0).drop_duplicates()
            ipgs.append(ipg)
            handle.close()

    # Concatenate all batches
    ipgs = pd.concat(ipgs).reset_index(drop=True).sort_values('id').drop_duplicates()
    found   = set(queries).intersection(set(ipgs.pid.unique()))
    missing = set(queries) - set(found)
    if verbose:
        print(f'Efetch: {len(ipgs)} IPG rows fetched, {len(found)} queries found, {len(missing)} missing accessions.', file=sys.stderr)

    # Filter error messages
    ipgs = ipgs[~ipgs.id.apply(lambda x: 'Cannot determine Ipg for accession' in str(x))]
    if ipgs.empty:
        ncbi.missing(missing)
        return pd.DataFrame(columns=[*cols, 'is_query', 'representative'])

    # Register query proteins
    ipgs['is_query'] = ipgs.pid.isin(queries).astype(int)

    # Process batches of different sizes
    if len(queries) == 1: # One query
        if not ipgs.empty:
            ipgs['representative'] = queries[0]

    # More than one query?
    else:
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
                ipg = ncbi.read('ipg', fetch=fetch, verbose=verbose)
                if ipgs.empty:
                    ipgs = ipg
                else:
                    ipgs = ipgs.append(ipg, ignore_index=True).drop_duplicates()
            ncbi.submit(queries) # Reset query list

    # Set id to numeric and update list of missing queries
    if not ipgs.empty:
        ipgs['id'] = pd.to_numeric(ipgs.id)
    missing = set([x for x in queries if x not in ipgs.pid.append(ipgs.representative).unique()])
    ncbi.missing(missing)
    if verbose:
        found = set([x for x in queries if x not in missing])
        print(f'{__name__}: {len(ipgs)} IPG rows fetched, {len(found)} queries found, {len(missing)} missing accessions.', file=sys.stderr)

    # Use assembly_reports to set priority
    if isinstance(assembly_reports,pd.DataFrame) or assembly_reports == True:
        col = ['wgs_master','ftp_path','isolate','paired_asm_comp','organism_name','submitter','loaded_from','infraspecific_name','asm_name','relation_to_type_material']
        if not isinstance(assembly_reports,pd.DataFrame):
            url = f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS'
            assembly_reports = type(ncbi)(query=list(ipgs[~ipgs.assembly.isna()].assembly.unique())).read('assembly_reports', baseurl=url, columns=col)
        ipgs = ipgs.merge(assembly_reports.drop(col, axis=1), left_on='assembly', right_on='assembly', how='left')

    # Return the expected dataframe
    return ipgs

# Load taxonomy data as a simple dataframe (tree as linearized path)
def taxonomy(ncbi, fetch=['ete3','entrez'], missing=False, ete3=None, preferred_taxa=None, verbose=False, *args, **kwargs):
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
      Columns: taxid, organism, domain, lineage, taxonomy

    Parameters:
      fetch    : list of fetchers
      missing  : retrieve missing entries only
      ete3     : ete3's NCBITaxa object
                 If set to true, a new NCBITaxa object is created
    '''

    # Set log format
    logger = logging.getLogger('rotifer.db.ncbi')
    if verbose:
        logger.setLevel(logging.DEBUG)
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
            tax = __taxonomy_from_ete3(ncbi, ete3=ete3, preferred_taxa=preferred_taxa, verbose=verbose)
        elif method == 'entrez':
            continue
        missing = ncbi.missing().copy()
        stack.append(tax)
        if verbose:
            print(f'{__name__}: {len(tax)} taxa loaded by {method}.', file=sys.stderr)
    if len(stack) > 0:
        tax = pd.concat(stack).drop_duplicates()
    else:
        tax = pd.DataFrame(columns='taxid organism domain lineage taxonomy'.split(' '))

    # Return
    ncbi.submit(queries)
    ncbi.missing(missing)
    if verbose:
        logger.info(f'main: {len(tax.taxid.unique())} taxids found and {len(missing)} missing!')
    return tax

# Internal methods

def __taxonomy_from_ete3(ncbi, ete3=None, preferred_taxa=None, verbose=False):
    # Make sure ete3 is a ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa
    if not isinstance(ete3,NCBITaxa):
        ete3 = NCBITaxa()
    if verbose:
        print(f'{__name__}: Searching for {len(ncbi.submit())} taxa.', file=sys.stderr)

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
    if verbose:
        print(f'{__name__}: loaded {len(s)} lineages.', file=sys.stderr)

    # Fetch names
    l = ete3.get_taxid_translator(list(l))
    if verbose:
        print(f'{__name__}: loaded {len(l)} taxon names.', file=sys.stderr)

    # Translate all lineages
    cols = ['taxid','organism','domain','taxonomy']
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
    if verbose:
        print(f'{__name__}: translated {len(data[cols[0]])} lineages.', file=sys.stderr)

    # Cleanup: remove missing taxids and register
    data = pd.DataFrame(data)
    from rotifer.taxonomy.utils import lineage
    data.insert(3, "lineage", lineage(data.taxonomy, preferred_taxa=preferred_taxa))
    ncbi.missing(missing)
    return data

# Is this library being used as a script?
if __name__ == '__main__':
    pass
