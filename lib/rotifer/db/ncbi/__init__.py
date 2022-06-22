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
import pandas as pd
import numpy as np

# Load NCBI configuration
import socket
from rotifer.core.functions import loadConfig as loadConfig
NcbiConfig = loadConfig(__name__.replace("rotifer.",":"))
if 'email' not in NcbiConfig:
    NcbiConfig['email'] = os.environ['USER'] + '@' + socket.gethostname()
if 'ftpserver' not in NcbiConfig:
    NcbiConfig['ftpserver'] = 'ftp.ncbi.nlm.nih.gov'
if 'NCBI_API_KEY' in os.environ:
    NcbiConfig['api_key'] = os.environ['NCBI_API_KEY']
else:
    NcbiConfig['api_key'] = None

# Import submodules
from rotifer.db.ncbi.ncbi import ncbi
from .ftp import *
from .entrez import *

# Controlling what can be exported with *
#
# This also ensure all documentation from imported
# submodules can be read here
#
# Note: to avoid exposing a function when * is used
# its name should start with a _
__all__ = [
        'assembly_reports',
        'neighbors',
        'ftp_get',
        'ftp_open',
        'elink',
        ]

# FUNCTIONS

# Gather neighbors of target genes by gene or gene product identifier
def neighbors(query=[], column='pid', assembly_reports=None, ipgs=None, exclude_type=['source','gene'], tries=1, verbose=0, *args, **kwargs):
    """
    Fetch gene neighborhoods directly from NCBI.

    Usage:
      # Simplest example:
      #   - single accession

      import pandas as pd
      import rotifer.db.ncbi as ncbi
      n = pd.concat(ncbi.neighbors(['WP_011017450.1']))

      # Slightly faster and more complicated:
      #   - two targets
      #   - pre-loaded assembly reports and IPGs
      #   - more neighbors (15)

      import pandas as pd
      import rotifer.db.ncbi as ncbi
      ar = ncbi.assemblies(taxonomy=True)
      i = ncbi.ncbi(["WP_010887045.1", "WP_011017450.1"]).read("ipg")
      n = pd.concat(
          ncbi.neighbors(
              ["WP_011017450.1","WP_010887045.1"],
              assembly_reports=ar,
              ipgs=i,
              before=15,
              after=15
          )
      )

    Returns:
      A generator of rotifer.genome.data.NeighborhoodDf dataframes

    Parameters:
      query  : list of strings
          List of accessions

      column : string
          Name of the column to scan for matches to the accessions
          See rotifer.genome.data.NeighborhoodDF

      assembly_reports : Pandas dataframe
        If not given, assembly reports are downloaded.
        See rotifer.db.ncbi.assemblies

      ipgs : rotifer.db.ncbi.read.ipg dataframe
        This parameter may be used to avoid downloading IPGs
        from NCBI. Example:

             import rotifer.db.ncbi as ncbiModule
             from rotifer.db.ncbi import ncbi
             i = ncbi(['WP_063732599.1']).read("ipg")
             n = ncbiModule.neighbors(['WP_063732599.1'], ipgs=i)

             Make sure it has the same colums as named
             by the ipg method in rotifer.db.ncbi.read

      exclude_type : list of strings
        Exclude rows by type (column 'type')

      tries : integer
        Number of attempts to download data

      verbose : integer
        Control verbosity for debugging

      Additional arguments are passed to the neighbors method
      of the rotifer.genome.data.NeighborhoodDF class
    """
    from Bio import SeqIO, Entrez
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
                assembly_reports = assemblies(verbose=verbose)
                attempt = tries + 1
            except:
                if verbose > 0:
                    msg = f'{__fn}: Failed to download assembly reports, {tries - attempt - 1} attempts left. Error: '
                    print(msg + str(sys.exc_info()[0]), file=sys.stderr)
                attempt += 1
        if not isinstance(assembly_reports,pd.DataFrame) or assembly_reports.empty:
            if verbose > 0:
                print(f'{__fn}: Failed to download assembly reports after {attempt} attempts.', file=sys.stderr)
            yield pd.DataFrame()

    # Fetch IPGs and discard those unrelated to the query
    if not isinstance(ipgs,pd.DataFrame):
        ipgs = ncbi(query).read('ipg', verbose=verbose)
        if ipgs.empty:
            print(f'{__fn}: failed to download IPGs for {len(query)} queries: {query}', file=sys.stderr)
            yield pd.DataFrame() # At his point, I can't handle missing IPGs for real NCBI protein accessions: fix using elink or efetch

    # Make sure we only process the IPGs of our queries
    selected = ipgs[ipgs.pid.isin(query) | ipgs.representative.isin(query)].id.unique()
    if not selected.any():
        if verbose > 0:
            print(f'{__fn}: No query was found in {len(ipgs.id.unique())} IPGs. Queries: {query}', file=sys.stderr)
        yield pd.DataFrame() # At his point, I can't handle missing IPGs for real NCBI protein accessions: fix using elink or efetch
    ipgs = ipgs[ipgs.id.isin(selected)]

    # Select best genomes per target: needs improvement!!!!!
    #
    # This code selects only one assembly per genome: identical
    # sequences will not be ignored!
    #
    # See Aureliano's code for better way (rotifer.db.ncbi.read.ipg)
    replaced = pd.Series(ipgs.representative.values, index=ipgs.pid).to_dict()
    in_ipg = set(ipgs.pid).union(set(ipgs.representative))
    found = set(query).intersection(in_ipg)
    missing = set(query) - found
    ipgs.sort_values(['id','order'], inplace=True)
    selected = ipgs.drop_duplicates(['id'], keep='first', ignore_index=True)
    if verbose > 0:
        print(f'{__fn}: {len(found)} queries were found in {len(selected.id.unique())} IPGs. Found: {found}, missing: {missing}', file=sys.stderr)

    # Download assemblies from NCBI's FTP site and process each of them
    failed = True
    pos = pd.Series(range(0,len(selected)))
    for s in pos:
        # Fetching next batch
        row = selected.iloc[s]
        acc = row.assembly
        if pd.isna(acc):
            acc = row.nucleotide
            if pd.isna(acc):
                continue # No solution yet for IPGs without nucleotides
            else:
                ndf = SeqIO.parse(Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=acc),"genbank")
        else:
            ndf = ncbi(acc).parse('genomes', assembly_reports=assembly_reports)

        # Parsing
        try:
            ndf = seqrecords_to_dataframe(ndf, exclude_type=exclude_type)
        except:
            if verbose > 0:
                msg = f'{__fn}: Failed to parse accession {acc}, {tries - pos.count(s)} attempts left. Error: '
                print(msg + sys.exc_info()[0], file=sys.stderr)
            if (pos.values == s).sum() < tries:
                pos.loc[pos.index.max() + 1] = s
            continue

        # Checking
        if ndf.empty:
            if verbose > 0:
                print(f'{__fn}: Empty NeighborhoodDF for accession {acc}, ignoring...', file=sys.stderr)
            continue
        elif column not in ndf.columns:
            if verbose > 0:
                print(f'{__fn}: NeighborhoodDF missing column {column} for accession {acc}', file=sys.stderr)
            continue

        # Identify queries
        select  = ndf[column].isin(in_ipg)
        if not select.any():
            if verbose > 0:
                print(f'{__fn}: No matches in {acc}. Ignoring accession...', file=sys.stderr)
            continue

        # Collecting neighbors
        #ndf = ndf.vicinity(select, *args, **kwargs)
        ndf = ndf.neighbors(select, *args, **kwargs)
        ndf['replaced'] = ndf.pid.replace(replaced)
        if 'block_id' in ndf.columns:
            kwargs['min_block_id'] = ndf.block_id.max() + 1
        else:
            kwargs['min_block_id'] = 1
        if not ndf.empty:
            failed = False
        yield ndf

    # I tried everything and found no queries in the target genomes
    if failed:
        print(f'{__fn}: No query was found after processing genomes in {len(ipgs.id.unique())} IPGs. Queries: {query}', file=sys.stderr)
        yield pd.DataFrame()

# Load NCBI assembly reports
def assemblies(baseurl=f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS', taxonomy=None, verbose=0):
    '''
    Load a table documenting all NCBI genome assemblies.

    By default, these concatenated tables are downloaded
    from the genomes/ASSEMBLY_REPORTS directory at NCBI's
    FTP site.

    Usage:
      # download from NCBI's FTP site
      import rotifer.db.ncbi as ncbi
      a = ncbi.assembly_reports()

      # Load local files at /db/ncbi
      b = ncbi.assembly_reports(baseurl="/db/ncbi")

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
    df = list()
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
        df.append(_)
        if verbose:
            logger.info(f'{__name__}: {url}, {len(_)} rows, {len(df)} loaded')
    df = pd.concat(df, ignore_index=True)
    if verbose:
        print(f'{__name__}: loaded {len(df)} assembly summaries.', file=sys.stderr)

    # Make sure the ftp_path columns refers to the ftp site as we expect
    if 'ftp_path' in df.columns:
        df.ftp_path = df.ftp_path.str.replace('https','ftp')

    # Add taxonomy
    if isinstance(taxonomy,pd.DataFrame) or taxonomy:
        if not isinstance(taxonomy,pd.DataFrame):
            from rotifer.db.ncbi import ncbi as ncbiClass
            ncbi = ncbiClass(df.taxid.unique().tolist())
            taxonomy = ncbi.read('taxonomy', ete3=taxonomy, verbose=verbose)
        if isinstance(taxonomy,pd.DataFrame):
            df = df.merge(taxonomy, left_on='taxid', right_on='taxid', how='left')
        if verbose:
            print(f'{__name__}: {len(df)} df left-merged with taxonomy dataframe.', file=sys.stderr)

    # Reset ncbi object, update missing list and return
    if verbose:
        logger.info(f'main: {len(df)} assembly reports loaded!')
    return df

# END
if __name__ == '__main__':
    pass

