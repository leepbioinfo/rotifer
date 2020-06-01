#!/usr/bin/env python3

import os
import sys
import logging
from rotifer.db.ncbi import NcbiConfig

# Load NCBI assembly reports
def assembly_reports(ncbi, baseurl=None, columns=[], query_type='assembly_accession', *args, **kwargs):
    '''
    Load NCBI assembly reports from a local directory or FTP.

    Usage:
      a = ncbi()
      b = a.read(method='assembly_reports')

      or, with filters,

      a = ncbi()
      a.submit(['GCF_900504695.1', 'GCF_004636045.1'])
      b = a.read(method='assembly_reports', columns=['id','accession'])

    Returns:
      Pandas DataFrame

    Parameters:
      baseurl    : URL or directory with assembly_summary_*.txt files
      columns    : (optional) list of columns to retrieve
      query_type : column name to use for filtering returned rows
                   Set to None to use pd.DataFrame.query

    Extra attributes:
      loaded_from : data source (same as baseurl)
    '''

    # Method dependencies
    import pandas as pd
    from glob import glob

    # Set log format
    logger = logging.getLogger('rotifer.db.ncbi')
    logger.setLevel(logging.DEBUG)
    logger.info(f'main: loading assembly reports...')

    # Choose URL
    if not baseurl:
        baseurl = f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS'
        if ('ROTIFER_DATA' in os.environ):
            localdir = os.path.join(os.environ["ROTIFER_DATA"],"genomes","ASSEMBLY_REPORTS")
            if os.path.exists(localdir) and glob(os.path.join(localdir,'assembly_summary_*.txt')):
                baseurl = localdir

    # Load assembly reports
    assemblies = list()
    for x in ['refseq', 'genbank', 'refseq_historical', 'genbank_historical']:
        if os.path.exists(baseurl): # Local file
            url = os.path.join(baseurl, f'assembly_summary_{x}.txt')
            if not os.path.exists(url):
                logger.warning(f'{__name__}: {url} not found. Ignoring...')
                continue
        else: # FTP
            url = f'{baseurl}/assembly_summary_{x}.txt'
        _ = pd.read_csv(url, sep ="\t", skiprows=[0])
        _.rename({'# assembly_accession':'assembly_accession'}, axis=1, inplace=True)
        _['source'] = x
        _['loaded_from'] = url
        assemblies.append(_)
        logger.info(f'{__name__}: {url}, {len(_)} rows, {len(assemblies)} loaded')
    assemblies = pd.concat(assemblies, ignore_index=True)

    # Filter rows and columns in the assemblies dataframe
    queries = ncbi.submit()
    if queries:
        if query_type in assemblies.columns: # Exact match values in a column
            assemblies = assemblies[assemblies[query_type].isin(queries)]
        else: # Query are logical statements
            for query in queries:
                assemblies = assemblies.query(query)
    ncbi.missing([ x for x in queries if x not in list(assemblies.assembly_accession) ])

    # Filter columns and return pandas object
    if columns:
        assemblies = assemblies.filter(columns)
    logger.info(f'main: {len(assemblies)} assembly reports loaded!')
    return assemblies

# Load Identical Protein Reports
def ipg(ncbi, fetch=['entrez'], assembly_reports=False, verbose=False, *args, **kwargs):
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
    cols = ['id','source','nucleotide','start','stop','strand','accession','description','organism','strain','assembly']
    ipgs = pd.DataFrame(columns=cols)

    # Set log format
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    logger.info(f'main: downloading IPG reports for {len(ncbi)} protein accessions...')

    # Backup and check queries
    queries = ncbi.submit()
    if len(queries) == 0:
        return pd.DataFrame(columns=[*cols, 'is_query', 'representative'])
    if verbose:
        print(f'{__name__}: searching {len(ncbi)} accessions...', file=sys.stderr)

    # Fetch using Epost+Efetch (recommended by NCBI but do we need it?)
    handle  = None
    missing = queries
    found   = []
    if len(queries) > 50:
        try:
            request = Entrez.epost(db = 'protein', api_key = ncbi.api_key(), id = ",".join(queries))
            result = Entrez.read(request)
            handle = Entrez.efetch(db='protein', rettype='ipg', retmode='text', api_key=ncbi.api_key(),
                                    webenv=result['WebEnv'], query_key=result['QueryKey'])
        except RuntimeError:
            if verbose:
                print(f'Epost+Efetch, {len(ncbi)} accessions. Runtime error: '+str(sys.exc_info()[1]), file=sys.stderr)
        except:
            if verbose:
                print(f'Epost+Efetch, {len(ncbi)} accessions. Unexpected error: '+str(sys.exc_info()[0:2]), file=sys.stderr)

        # Parse fetched IPG reports
        if handle:
            ipgs = pd.read_csv(handle, sep='\t', names=cols, header=0).drop_duplicates()
            found   = set([x for x in queries if x in ipgs['accession'].unique()])
            missing = set([x for x in queries if x not in found])
            handle.close() # Close Bio.Entrez _io.TextWrapper handle

        # Report on progress
        if verbose:
            print(f'Epost+Efetch, {len(ipgs)} rows in IPG report, {len(found)} queries found, {len(missing)} missing accessions', file=sys.stderr)

    # Try using Efetch directly!
    if missing:
        handle = None
        try:
            handle = Entrez.efetch(db='protein', rettype='ipg', retmode='text', api_key=ncbi.api_key(), id = ",".join(missing))
        except RuntimeError:
            if verbose:
                print(f'Efetch: {len(queries)} queries, {len(missing)} accessions. Runtime error: '+str(sys.exc_info()[1]), file=sys.stderr)
        except:
            if verbose:
                print(f'Efetch: {len(queries)} queries, {len(missing)} accessions. Unexpected error: '+str(sys.exc_info()[0:2]), file=sys.stderr)

        # Merge Efetch results with Epost+Efetch results
        if handle:
            ipg = pd.read_csv(handle, sep='\t', names=cols, header=0).drop_duplicates()
            if ipgs.empty:
                ipgs = ipg
            else:
                ipgs = ipgs.append(ipg, ignore_index=True).drop_duplicates()
            found   = set([x for x in queries if x in ipgs['accession'].unique()])
            missing = set([x for x in queries if x not in found])
            handle.close()

        # Report progress
        if verbose:
            print(f'Efetch: {len(ipgs)} IPG rows fetched, {len(found)} queries found, {len(missing)} missing accessions.', file=sys.stderr)

    # Filter error messages
    ipgs = ipgs[~ipgs.id.apply(lambda x: 'Cannot determine Ipg for accession' in str(x))]
    if ipgs.empty:
        ncbi.missing(missing)
        return pd.DataFrame(columns=[*cols, 'is_query', 'representative'])

    # Register query proteins
    ipgs['is_query'] = ipgs.accession.isin(queries).astype(int)

    # Process batches of different sizes
    if len(queries) == 1: # One query
        if not ipgs.empty:
            ipgs['representative'] = queries[0]

    # More than one query?
    else:
        #  Register first query protein as representative
        rep = ipgs.loc[ipgs['is_query'] == 1,['id','accession']].drop_duplicates('id', keep='first')
        idmap = {k:v for k,v in zip(rep['id'].values, rep['accession'].values) }
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
    missing = set([x for x in queries if x not in ipgs.accession.append(ipgs.representative).unique()])
    ncbi.missing(missing)
    if verbose:
        found = set([x for x in queries if x not in missing])
        print(f'{__name__}: {len(ipgs)} IPG rows fetched, {len(found)} queries found, {len(missing)} missing accessions.', file=sys.stderr)

    # Use assembly_reports to set priority
    if isinstance(assembly_reports,pd.DataFrame) or assembly_reports == True:
        col = ['assembly_accession','bioproject','biosample','refseq_category','taxid','species_taxid',
                'version_status','assembly_level','release_type','genome_rep','seq_rel_date','source']
        if not isinstance(assembly_reports,pd.DataFrame):
            url = f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS'
            assembly_reports = type(ncbi)(query=list(ipgs[~ipgs.assembly.isna()].assembly.unique())).read('assembly_reports', baseurl=url, columns=col)
        ipgs = ipgs.merge(assembly_reports[col], left_on='assembly', right_on='assembly_accession', how='left')

    return ipgs

# Is this library being used as a script?
if __name__ == '__main__':
    pass
