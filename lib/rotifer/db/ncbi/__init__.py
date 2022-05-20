# Copyright 2020 by Robson F. de Souza.  All rights reserved.
# This file is part of the Rotifer distribution and governed by 
# the "BSD 3-Clause License".
#
# Please see the LICENSE file that should have been included as part of this
# package.

r"""
=============================
Rotifer's NCBI database tools
=============================

Class attributes
----------------

NcbiConfig : NCBI configuration
  Automatically loaded from ~/.rotifer/etc/db/ncbi.yml
"""

# Import external modules
import os
import sys
import logging
import socket
#import types
import pandas as pd

# Load module's configuration and defaults
from rotifer.core import GlobalConfig
from rotifer.core.functions import loadConfig
#import rotifer.db.ncbi.ncbi as ncbi

# Load NCBI configuration
NcbiConfig = loadConfig(__name__.replace("rotifer.",":"))
if 'email' not in NcbiConfig:
    NcbiConfig['email'] = os.environ['USER'] + '@' + socket.gethostname()
if 'ftpserver' not in NcbiConfig:
    NcbiConfig['ftpserver'] = 'ftp.ncbi.nlm.nih.gov'
if 'NCBI_API_KEY' in os.environ:
    NcbiConfig['api_key'] = os.environ['NCBI_API_KEY']

# Import submodules
from rotifer.db.ncbi.ncbi import ncbi

# FUNCTIONS

# Gather neighbors of target genes by gene or gene product identifier
def neighbors(query=[], column='pid', assembly_reports=None, ipgs=None, exclude_type=['source','gene'], batch_size=1, tries=1, verbose=0, *args, **kwargs):
    """
    Fetch gene neighborhoods directly from NCBI.

    Usage:
      import rotifer.db.ncbi as ncbiClass
      n = ncbiClass.neighbors(['WP_063732599.1','WP_063732345.1'])

      # Using previously loaded assembly reports (slightly faster)
      import rotifer.db.ncbi as ncbiClass
      a = ncbiClass.assembly_reports()
      b = ncbiClass.neighbors(['WP_063732599.1'], assembly_reports=a)

    Returns:
      A rotifer.genome.data.NeighborhoodDf dataframe

    Parameters:
      query  : list of accessions
      column : column to scan for matches to accessions

      assembly_reports : rotifer.db.ncbi.assembly_reports dataframe
                         If not given, assembly reports are downloaded

                         Keep in mind that ONLY assemblies listed in this
                         dataframe are searched for matches to queries.

                         Therefore, one can use a slice to fetch entries from
                         a selected set of genomes:

                         import rotifer.db.ncbi as ncbiClass
                         a = ncbiClass.assembly_reports()
                         a = a[a.assembly.isin(['GCF_001650215.1'])
                         b = ncbiClass.neighbors(['WP_063732599.1'], assembly_reports=a)

      ipgs : rotifer.db.ncbi.read.ipg dataframe
             This parameter may be used to avoid downloading IPGs
             from NCBI. Example:

             import rotifer.db.ncbi as ncbiClass
             from rotifer.db.ncbi import ncbi
             i = ncbi(['WP_063732599.1']).read("ipg")
             n = ncbiClass.neighbors(['WP_063732599.1'], ipgs=i)

             Make sure it has the same colums as named
             by the ipg method in rotifer.db.ncbi.read

      exclude_type : exclude features by type

      batch_size : how many genomes to download and process at each
                   from the NCBI FTP site at each parsing round

      tries : number of attempts to download data

      verbose : (integer) control verbosity for debugging

      Other arguments are supported by neighbors from rotifer.genome.data
    """

    import numpy as np
    from rotifer.db.ncbi import ncbi
    from rotifer.genome.utils import seqrecords_to_dataframe
    __fn = "rotifer.db.ncbi.neighbors"

    # Adjust minimum block ID
    if not "min_block_id" in kwargs:
        kwargs["min_block_id"] = 1

    # Prepare a rotifer.db.ncbi object
    if not isinstance(query,list):
        if isinstance(query, pd.Series):
            query = query.tolist()
        else:
            query = [ query ]

    # Fetch assembly reports
    if not isinstance(assembly_reports,pd.DataFrame) or assembly_reports.empty:
        attempt = 0
        while attempt < tries:
            try:
                assembly_reports = assembly_reports(verbose=verbose)
                attempt = tries + 1
            except:
                if verbose > 0:
                    msg = f'{__fn}: Failed to download assembly reports, {tries - attempt - 1} attempts left. Error: '
                    print(msg + str(sys.exc_info()[0]), file=sys.stderr)
                attempt += 1
        if not isinstance(assembly_reports,pd.DataFrame) or assembly_reports.empty:
                if verbose > 0:
                    print(f'{__fn}: Failed to download assembly reports after {attempt + 1} attempts left.', file=sys.stderr)
                return pd.DataFrame()

    # Fetch IPGs
    if not isinstance(ipgs,pd.DataFrame):
        ncbiObj = ncbi(query)
        try:
            ipgs = ncbiObj.read('ipg', verbose=verbose)
        except:
            if verbose > 0:
                print(f'{__fn}: unexpected error while downloading IPGs: '+str(sys.exc_info()[0:2]), file=sys.stderr)
                return None

    # Filter IPGs based on assembly reports
    ipgs = ipgs[ipgs.assembly.isin(assembly_reports.assembly)]

    # IPGs are valid if pid or representative is a query
    found = ipgs.pid.isin(query)
    found = found | ((~ipgs.id.isin(ipgs[found].id)) & ipgs.representative.isin(query))
    ipgs  = ipgs[found]
    if ipgs.empty:
        if verbose > 0:
            print(f'{__fn}: Empty IPG dataframe, no IPG reports found! Aborting...', file=sys.stderr)
            return None

    # Select best genomes per target: needs improvement!!!!!
    #
    # This code selects only one assembly per genome: identical
    # sequences will not be ignored!
    #
    # See Aureliano's code for better way (rotifer.db.ncbi.read.ipg)
    ipgs = ipgs.groupby(['id']).agg({'pid':'first','assembly':'first','representative':'first'}).reset_index()
    found = set(ipgs.pid)
    missing = set(query) - found - set(ipgs.representative)
    if verbose > 0:
        print(f'{__fn}: {len(found)} queries were found in {len(ipgs)} IPGs. Found: {found}, missing: {missing}', file=sys.stderr)

    # Prepare list of assemblies (use first found in IPG)
    # Download and parse assemblies
    ndf = []
    assemblies = ipgs.assembly.sort_values().unique().tolist()
    pos = list(range(0,len(assemblies),batch_size))
    for s in pos:
        # Fetching next batch
        e = s + batch_size
        ids = assemblies[s:e]
        genomes = ncbi(ids).parse('genomes', assembly_reports=assembly_reports)

        # Parsing
        try:
            genomes = seqrecords_to_dataframe(genomes, exclude_type=exclude_type)
        except:
            if verbose > 0:
                print(f'{__fn}: Failed to parse batch {s}, {tries - pos.count(s)} attempts left. Genomes: {ids}. Error: '+sys.exc_info()[0], file=sys.stderr)
            if pos.count(s) < tries:
                pos.append(s)
            continue

        # Checking
        if genomes.empty:
            if verbose > 0:
                print(f'{__fn}: Empty NeighborhoodDF for batch {s}, ignoring assemblies {ids}', file=sys.stderr)
            continue
        elif column not in genomes.columns:
            if verbose > 0:
                print(f'{__fn}: NeighborhoodDF missing column {column} at batch {s}, assemblies {ids}', file=sys.stderr)
            continue

        # Searching targets
        select  = genomes[column].isin(found)
        if not select.any():
            if verbose > 0:
                print(f'{__fn}: No matches in {", ".join(genomes.assembly.unique())}. Ignoring batch {s}:{e} ...', file=sys.stderr)
                continue

        # Collecting neighbors
        genomes = genomes.neighbors(select, *args, **kwargs)
        #genomes = genomes.vicinity(select, *args, **kwargs)
        ndf.append(genomes)
        kwargs["min_block_id"] = genomes.block_id.max() + 1

    # Merging neighborhoods
    if len(ndf) == 0: return pd.DataFrame()
    ndf = pd.concat(ndf)
    ndf.drop_duplicates(inplace=True)
    ndf.reset_index(drop=True, inplace=True)
    replaced = pd.Series(ipgs.representative.values, index=ipgs.pid).to_dict()
    ndf['replaced'] = ndf.pid.replace(replaced)

    return ndf

# Load NCBI assembly reports
def assembly_reports(baseurl=f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS', taxonomy=None, verbose=0):
    '''
    Load NCBI assembly reports from a local directory or from the NCBI FTP site.

    Usage:
      # download from NCBI's FTP site
      from rotifer.db.ncbi import assembly_reports
      a = assembly_reports()

      # Load local files at /db/ncbi
      b = assembly_reports(baseurl="/db/ncbi")

    Returns:
      Pandas DataFrame

    Parameters:
      baseurl    : URL or directory with assembly_summary_*.txt files
      taxonomy   : ete3's NCBITaxa object
                   If set to true, a new NCBITaxa object is created

    Extra columns added by this method:
      source : NCBI's source database
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

    # Add taxonomy
    if isinstance(taxonomy,pd.DataFrame) or taxonomy:
        if not isinstance(taxonomy,pd.DataFrame):
            from rotifer.db.ncbi import ncbi as ncbiClass
            ncbi = ncbiClass(assemblies.taxid.unique().tolist())
            taxonomy = ncbi.read('taxonomy', ete3=taxonomy, verbose=verbose)
        if isinstance(taxonomy,pd.DataFrame):
            assemblies = assemblies.merge(taxonomy, left_on='taxid', right_on='taxid', how='left')
        if verbose:
            print(f'{__name__}: {len(assemblies)} assemblies left-merged with taxonomy dataframe.', file=sys.stderr)

    # Filter columns and return pandas object
    if columns:
        assemblies = assemblies.filter(columns)

    # Reset ncbi object, update missing list and return
    if verbose:
        logger.info(f'main: {len(assemblies)} assembly reports loaded!')
    return assemblies

# END
if __name__ == '__main__':
    pass

