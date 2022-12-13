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

config : NCBI configuration
  Automatically loaded from ~/.rotifer/etc/db/ncbi.yml
"""

# Import external modules
import os
import sys
import types
import socket
import typing
import logging
import importlib
import numpy as np
import pandas as pd
from tqdm import tqdm
from copy import deepcopy

# Import rotifer modules
import rotifer
from rotifer import GlobalConfig
logger = rotifer.logging.getLogger(__name__)

# Load NCBI configuration
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.delegator
from rotifer.core.functions import loadConfig
config = loadConfig(
    __name__.replace("rotifer.",":"),
    defaults = {
        "mirror": os.path.join(os.environ["ROTIFER_DATA"] if 'ROTIFER_DATA' in os.environ else "/databases","genomes"),
        'email': os.environ['USER'] + '@' + socket.gethostname(),
        'ftpserver': 'ftp.ncbi.nlm.nih.gov',
        'api_key': os.environ['NCBI_API_KEY'] if 'NCBI_API_KEY' in os.environ else None,
        'cursor_methods': {
            'entrez': 'rotifer.db.ncbi.entrez',
            'ete3': 'rotifer.db.local.ete3',
            'ftp': 'rotifer.db.ncbi.ftp',
            'mirror': 'rotifer.db.ncbi.mirror',
            'sqlite3': 'rotifer.db.sql.sqlite3',
        }
    }
)
NcbiConfig = config # for compatibility but deprecated: to be removed!

# Classes

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    def __init__(
            self,
            methods=['mirror','ftp'],
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=10,
            timeout=10,
            basepath = config["mirror"],
            cache=GlobalConfig['cache'],
            *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads','cache','basepath']
        self.sleep_between_tries = sleep_between_tries
        self.timeout = timeout
        self.basepath = basepath
        self.cache = cache
        super().__init__(methods=methods, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

class GenomeFeaturesCursor(rotifer.db.methods.GenomeFeaturesCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes

    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> from rotifer.db.ncbi as ncbi
    >>> gfc = ncbi.GenomeFeaturesCursor(progress=True)
    >>> df = gfc.fetchall(g)

    Parameters
    ----------
    exclude_type: list of strings
      List of names for the features that must be ignored
    autopid: boolean
      Automatically set protein identifiers
    codontable: string por int, default 'Bacterial'
      Default codon table, if not set within the data
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    threads: integer, default 15
      Number of processes to run parallel downloads
    batch_size: int, default 1
      Number of accessions per batch
    cache: path-like string
      Where to place temporary files

    """
    def __init__(
            self,
            methods=['mirror','ftp'],
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            sleep_between_tries=1,
            batch_size=None,
            threads=15,
            timeout=10,
            basepath = config["mirror"],
            cache=GlobalConfig['cache'],
            *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads','cache','basepath']
        self.sleep_between_tries = sleep_between_tries
        self.timeout = timeout
        self.basepath = basepath
        self.cache = cache
        super().__init__(methods=methods, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

class GeneNeighborhoodCursor(rotifer.db.core.BaseGeneNeighborhoodCursor):
    """
    Fetch gene neighborhoods as dataframes
    ======================================

    This class implements the logic to search for genomic
    patches centered around a target coding gene.

    The target genes must be identified based on the accession
    numbers of their proteins.

    For the multi-query methods (fetchone and fetchall)
    results always return in a random order.

    Usage
    -----
    Using the dictionary-like interface, fetch the gene
    neighborhood around the gene encoding a target protein:

    >>> import rotifer.db.ncbi a ncbi
    >>> gnc = ncbi.GeneNeighborhoodCursor()
    >>> df = gnc["EEE9598493.1"]

    Fetch all gene neighborhoods for a sample of proteins:

    >>> q = ['WP_012291365.1','WP_013208129.1','WP_122330970.1']
    >>> df = gnc.fetchall(q)

    Process gene neighborhoods while downloading:

    >>> for n in gnc.fetchone(q):
    >>>     do_something(n)

    Parameters
    ----------
    *Important note:*

    All parameters for initialization of this class are acessible
    as mutable attributes and can be modified to tune the cursor's
    behaviour.

    column : string
      Name of the column to scan for matches to the accessions
      See rotifer.genome.data.NeighborhoodDF
    before : int
      Keep at most this number of features, of the same type as the
      target, before each target
    after  : nt
      Keep at most this number of features, of the same type as the
      target, after each target
    min_block_distance : int
      Minimum distance between two consecutive blocks
    strand : string
      How to evaluate rows concerning the value of the strand column
      Possible values for this option are:

      - None : ignore strand
      - same : same strand as the targets
      -    + : positive strand features and targets only
      -    - : negative strand features and targets only

    fttype : string
      How to process feature types of neighbors
      Supported values:

      - same : consider only features of the same type as the target
      - any  : ignore feature type and count all features when
               setting neighborhood boundaries

    eukaryotes : boolean, default False
      If set to True, neighborhood data for eukaryotic genomes
    save : string, default None
      If set, save processed batches to the path given
    replace : boolean, default True
      When save is set, whether to replace that file
    mirror : path
      Path to a mirror of the NCBI's FTP genomes directory.
    exclude_type: list of strings
      List of names for the features that must be ignored
    autopid: boolean
      Automatically set protein identifiers
    codontable: string por int, default 'Bacterial'
      Default codon table, if not set within the data
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    threads: integer, default 15
      Number of processes to run parallel downloads
    batch_size: int, default 1
      Number of accessions per batch
    cache: path-like string
      Where to place temporary files

    Internal state attributes
    -------------------------
    Objects of this class modify the folllowing read/write
    attributes when fetch methods are called:
    * missing
      A Pandas DataFrame describing errors and messages for
      failed attempts to download gene neighborhoods.

    """
    def __init__(
            self,
            column = 'pid',
            before = 7,
            after = 7,
            min_block_distance = 0,
            strand = None,
            fttype = 'same',
            eukaryotes=False,
            save=None,
            replace=False,
            mirror=None,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            batch_size=None,
            threads=15,
            cache=GlobalConfig['cache'],
            *args, **kwargs
        ):
        from rotifer.db.ncbi import ftp
        from rotifer.db.ncbi import entrez

        # Setup special attributes
        self._shared_attributes = [
            'column','before','after','min_block_distance','strand','fttype','eukaryotes',
            'exclude_type','autopid','codontable',
            'progress','tries','batch_size','threads','cache',
        ]
        self.cursors = [
            ftp.GeneNeighborhoodCursor(),
            entrez.GeneNeighborhoodCursor()
        ]
        if mirror:
            from rotifer.db.ncbi import mirror as rdnm
            cursor = rdnm.GeneNeighborhoodCursor(basepath=mirror)
            self.cursors.insert(0,cursor)
        if save:
            from rotifer.db.sql import sqlite3 as rdss
            cursor = rdss.GeneNeighborhoodCursor(save, replace=replace)
            self.cursors.insert(0,cursor)

        # Setup simple attributes
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
        self.save = save
        self.eukaryotes = eukaryotes
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable
        self.progress = progress
        self.tries = tries
        self.batch_size = batch_size
        self.threads = threads
        self.cache = cache
        self.missing = pd.DataFrame(columns=["noipgs","eukaryote","assembly","error",'class'])

    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        if hasattr(self,'_shared_attributes') and name in self._shared_attributes:
            for cursor in self.cursors:
                if hasattr(cursor,name):
                    cursor.__setattr__(name,value)

    def __getitem__(self, protein, ipgs=None):
        """
        Dictionary-like access to gene neighbors.

        Usage
        -----
        >>> import rotifer.db.ncbi as ncbi
        >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
        >>> n = gnc["WP_063732599.1"]

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db.ncbi import entrez
          >>> import rotifer.db.ncbi as ncbi
          >>> ic = entrez.IPGCursor(batch_size=1)
          >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.__getitem__(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        Generator of rotifer.genome.data.NeighborhoodDF
         """
        from rotifer.genome.utils import seqrecords_to_dataframe
        result = seqrecords_to_dataframe([])
        for cursor in self.cursors:
            result = cursor.__getitem__(protein, ipgs=ipgs)
            if not isinstance(result,types.NoneType) and len(result) > 0:
                if self.save and (cursor != self.cursors[0]):
                    self.cursors[0].insert(result)
                break
        return result

    def fetchone(self, accessions, ipgs=None):
        """
        Fetch each gene neighborhood iteratively.

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db.ncbi import entrez
          >>> import rotifer.db.ncbi as ncbi
          >>> ic = entrez.IPGCursor(batch_size=1)
          >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> for x in gnc.fetchone(['WP_063732599.1'], ipgs=i):
          >>>    do_something(x)

        Returns
        -------
        Generator of rotifer.genome.data.NeighborhoodDF
        """
        from rotifer.genome.utils import seqrecords_to_dataframe

        # Copy identifiers and remove redundancy
        targets = deepcopy(accessions)
        if not isinstance(targets,typing.Iterable) or isinstance(targets,str):
            targets = [targets]
        targets = set(targets)
        todo = deepcopy(targets)

        # Make sure we have IPGs
        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            if self.progress:
                logger.warn(f'Downloading IPGs for {len(todo)} proteins...')
            size = self.batch_size
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries, threads=self.threads)
            ipgs = ic.fetchall(list(todo))
            if len(ic.missing):
                self.update_missing(ic.missing, np.nan, "No IPGs")
                todo -= ic.missing

        # Select IPGs corresponding to our queries
        ipgs = ipgs[ipgs.id.isin(ipgs[ipgs.pid.isin(todo) | ipgs.representative.isin(todo)].id)]
        missing = todo - set(ipgs.pid).union(ipgs.representative)
        if missing:
            self.update_missing(missing,np.NaN,"Not found in IPGs")
            todo = todo - missing
        if len(ipgs) == 0:
            return [seqrecords_to_dataframe([])]

        # Call cursors
        lost = 'noipgs == False'
        if not self.eukaryotes:
            lost += ' and eukaryote == False'
        for i in range(0,len(self.cursors)):
            cursor = self.cursors[i]
            if len(todo) == 0:
                break
            for result in cursor.fetchone(todo, ipgs=ipgs):
                if self.save and i > 0:
                    self.cursors[0].insert(result)
                found = self.getids(result, ipgs=ipgs)
                for c in [self] + self.cursors[0:i+1]:
                    c.missing.drop(found, axis=0, inplace=True, errors="ignore")
                for s in cursor.missing.iterrows():
                    if s[0] in targets:
                        self.missing.loc[s[0]] = s[1]
                todo = set(self.missing.query(lost).index)
                yield result

    def fetchall(self, proteins, ipgs=None):
        """
        Fetch all gene neighborhoods in a single dataframe.

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db.ncbi import entrez
          >>> import rotifer.db.ncbi as ncbi
          >>> ic = entrez.IPGCursor(batch_size=1)
          >>> gnc = ncbi.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.fetchall(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        from rotifer.genome.utils import seqrecords_to_dataframe
        stack = []
        for df in self.fetchone(proteins, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

class TaxonomyCursor(rotifer.db.delegator.SequentialDelegatorCursor):
    def __init__(self, methods=['ete3','entrez'], progress=True, tries=3, sleep_between_tries=1, batch_size=None, threads=10, *args, **kwargs):
        self._shared_attributes = ['progress','tries','sleep_between_tries','batch_size','threads']
        self.sleep_between_tries = sleep_between_tries
        super().__init__(methods=methods, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.taxcols = ['taxid','organism','superkingdom','lineage','classification','alternative_taxids']

    def getids(self, obj, *args, **kwargs):
        if not (isinstance(obj,list) or isinstance(obj,tuple)):
            obj = [ obj ]
        ids = set()
        for item in obj:
            ids.update(set(item.taxid.astype(str)))
            if 'alternative_taxids' in item:
                aids = item.alternative_taxids.dropna().astype(str)
                aids = aids.str.split(",").explode().dropna()
                ids.update(aids)
        return ids

    def __getitem__(self, accessions, *args, **kwargs):
        """
        Dictionary-like access to data.

        Usage
        -----
        >>> import rotifer.db.ncbi as ncbi
        >>> tc = ncbi.TaxonomyCursor(progress=True)
        >>> t = tc[2599]

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.

        Returns
        -------
        Pandas dataframe
        """
        result = super().__getitem__(accessions, *args, **kwargs)
        if len(result) == 0:
            return pd.DataFrame(columns=self.taxcols)
        elif isinstance(result,list):
            return pd.concat(result, ignore_index=True)
        else:
            return result

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch data for all accessions.

        Parameters
        ----------
        accessions: list of database identifiers
          Database identifiers.

        Returns
        -------
        Pandas dataframe
        """
        df = super().fetchall(accessions, *args, **kwargs)
        if len(df) == 0:
            return pd.DataFrame(columns=self.taxcols)
        else:
            return pd.concat(df, ignore_index=True)

# FUNCTIONS

# Load NCBI assembly reports
def assemblies(baseurl=f'ftp://{config["ftpserver"]}/genomes/ASSEMBLY_REPORTS', taxonomy=True, progress=True):
    '''
    Load a table documenting all NCBI genome assemblies.

    By default, these concatenated tables are downloaded
    from the genomes/ASSEMBLY_REPORTS directory at NCBI's
    FTP site.

    Usage
    -----
    Download from NCBI's FTP site

    >>> import rotifer.db.ncbi as ncbi
    >>> a = ncbi.assembly_reports()

    Load local files at /db/ncbi

    >>> b = ncbi.assembly_reports(baseurl="/db/ncbi")
      
    If working at NIH servers use

    >>> a = ncbi.assembly_reports(baseurl="/am/ftp-genomes/ASSEMBLY_REPORTS")

    Parameters
    ----------
    baseurl: string
      URL or directory with assembly_summary_*.txt files
    taxonomy: boolean, default True
      If set to true, taxonomy data is added to the table

    Returns
    -------
    Pandas DataFrame

    Notes
    -----
    Some columns are added to the original table:
    source : NCBI's source database
    loaded_from : data source (same as baseurl)
    '''

    # Method dependencies
    import pandas as pd
    from glob import glob
    origLevel = logger.getEffectiveLevel()
    if progress:
        rotifer.logger.setLevel(rotifer.logging.INFO)
    logger.info(f'main: loading assembly reports...')

    # Load assembly reports
    df = list()
    for x in ['refseq', 'genbank', 'refseq_historical', 'genbank_historical']:
        if os.path.exists(baseurl): # Local file
            url = os.path.join(baseurl, f'assembly_summary_{x}.txt')
            if not os.path.exists(url):
                logger.warning(f'{__name__}: {url} not found. Ignoring...')
                continue
        else: # FTP
            url = f'{baseurl}/assembly_summary_{x}.txt'
        _ = pd.read_csv(url, sep ="\t", skiprows=[0])
        _.rename({'# assembly_accession':'assembly'}, axis=1, inplace=True)
        _['source'] = x
        _['loaded_from'] = url
        df.append(_)
        logger.info(f'{url}, {len(_)} rows, {len(df)} loaded')
    df = pd.concat(df, ignore_index=True)
    df.taxid = df.taxid.astype(str)
    logger.info(f'loaded {len(df)} assembly summaries.')

    # Make sure the ftp_path columns refers to the ftp site as we expect
    if 'ftp_path' in df.columns:
        df.ftp_path = df.ftp_path.str.replace('https','ftp')

    # Add taxonomy
    if taxonomy:
        cursor = TaxonomyCursor(progress=progress)
        taxonomy = cursor.fetchall(df.taxid.unique().tolist())
        taxonomy['_same'] = (taxonomy.taxid == taxonomy.alternative_taxids).astype(int)
        taxonomy.sort_values(['taxid','_same'], ascending=True, inplace=True)
        taxonomy.drop('_same', axis=1, inplace=True)
        taxonomy.drop_duplicates('taxid', keep='first', inplace=True)
        df = df.merge(taxonomy, left_on='taxid', right_on='taxid', how='left')
        logger.info(f'{len(df)} df left-merged with taxonomy dataframe.')

    # Reset ncbi object, update missing list and return
    logger.info(f'main: {len(df)} assembly reports loaded!')
    if progress:
        rotifer.logger.setLevel(origLevel)
    return df

# END
if __name__ == '__main__':
    pass
