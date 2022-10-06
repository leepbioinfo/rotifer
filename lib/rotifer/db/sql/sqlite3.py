
__doc__ = """
Rotifer connections to SQL databases
====================================
"""

import re
import os
import types
import typing
import sqlite3
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import StringIO

import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
from rotifer.db.core import BaseCursor
from rotifer.genome.data import NeighborhoodDF
from rotifer.genome.utils import seqrecords_to_dataframe
import rotifer.devel.beta.sequence as rdbs
logger = rotifer.logging.getLogger(__name__)

# Defaults
_config = loadConfig(__name__.replace('rotifer.',':'))
if not _config:
    _config = {}
_config = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
    **_config
}

class GeneNeighborhoodCursor(BaseCursor):
    """
    Fetch gene neighborhoods from SQLite3 database
    ==============================================

    This class implements storage of fully defined
    gene neighborhoods in a local SQLite3 database

    Input data should be gene neighborhoods returned
    by remote providers, such as NCBI.

    Usage
    -----
    Using the dictionary-like interface, fetch the gene
    neighborhood around the gene encoding a target protein:

    >>> from rotifer.db.sql import sqlite3 as rdss
    >>> gnc = rdss.GeneNeighborhoodCursor("neighbors.sqlite3")
    >>> df = gnc["EEE9598493.1"]

    Fetch all gene neighborhoods for a sample of proteins:

    >>> q = ['WP_012291365.1','WP_013208129.1','WP_122330970.1']
    >>> df = gnc.fetchall(q)

    Process gene neighborhoods while loading:

    >>> for n in gnc.fetchone(q):
    >>>     do_something(n)

    Parameters
    ----------
    *Important note:*

    All parameters for initialization of this class are acessible
    as mutable attributes and can be modified to tune the cursor's
    behaviour.

    path: string
      Path to the localstorage
    format : string, default ```tsv```
      Format of the local storage.
      Should match one of following the supported backends:
      * tsv: TAB-separated text tables
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
    replace : boolean, default False
      When save is set, whether to replace that file
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
    Objects of this class modify two main read/write attributes
    when fetch methods are called:
    * missing
      A Pandas DataFrame describing errors and messages for
      failed attempts to download gene neighborhoods.

    """
    def __init__(
            self,
            path,
            replace = False,
            column = 'pid',
            before = 7,
            after = 7,
            min_block_distance = 0,
            strand = None,
            fttype = 'same',
            eukaryotes=False,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
            cache=GlobalConfig['cache'],
        ):
        self.path = path
        self.replace = replace
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
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
        if os.path.exists(self.path) and self.replace:
            os.remove(self.path)
        self._dbconn = sqlite3.connect(self.path)
        self._dbconn.execute('CREATE TEMPORARY TABLE IF NOT EXISTS queries (id)')

    def __setattr__(self,name,value):
        super().__setattr__(name,value)
        if name == "column":
            if isinstance(value,str) or not isinstance(value, typing.Iterable):
                super().__setattr__(name,[value])
            if "pid" in self.column and "replaced" not in self.column:
                self.column.append("replaced")
            if "replaced" in self.column and "pid" not in self.column:
                self.column.append("pid")

    def _add_to_missing(self, accessions, assembly, error):
        err = [False,False,assembly,error,__name__]
        if "Eukaryotic" in error:
            err[1] = True
        if "IPG" in error:
            err[0] = True
        if not isinstance(accessions,typing.Iterable) or isinstance(accessions,str):
            accessions = [accessions]
        for x in accessions:
            self.missing.loc[x] = err

    def _has_table(self, name):
        sql = self._dbconn.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{name}'").fetchall()
        return len(sql) > 0

    def _fetch_from_sql(self, accessions):
        if not self._has_table("neighborhoods"):
            return seqrecords_to_dataframe([])
        cursor = self._dbconn.cursor()
        cursor.executemany("INSERT INTO queries VALUES (?)",[ (x,) for x in list(accessions) ])
        self._dbconn.commit()
        sql = []
        for col in self.column:
            sql.append(f'SELECT n2.* FROM queries as q inner join neighborhoods as n on (q.id = n.{col}) inner join neighborhoods as n2 using (block_id)')
        sql = " UNION ".join(sql)
        df = NeighborhoodDF(pd.read_sql(sql, self._dbconn))
        self._dbconn.execute(f'DELETE FROM queries')
        return df

    def __getitem__(self, accession, ipgs=None):
        """
        Dictionary-like access to gene neighbors.
        """
        if not isinstance(accession, typing.Iterable) or isinstance(accession,str):
            accession = [accession]
        accession = set(accession)

        # Load stored data
        stored = self._fetch_from_sql(accession)
        if len(stored) == 0:
            self._add_to_missing(accession,np.NaN,"Not found in storage")
            return seqrecords_to_dataframe([])
        found = stored.melt(id_vars=["assembly"], value_vars=self.column, var_name="type", value_name="id")
        missing = accession - set(found.id)

        # Find
        if len(missing) > 0 and 'pid' in self.column and (not isinstance(ipgs,types.NoneType)) and len(ipgs) > 0:
            ipgs = ipgs.melt(id_vars=['id','assembly'], value_vars=['pid','representative'], var_name="type", value_name='pid')
            ipgs.drop_duplicates(inplace=True)
            ipgs.rename({'id':'ipg','pid':'id'}, axis=1, inplace=True)
            storedIPGs = ipgs[ipgs.ipg.isin(ipgs[ipgs.id.isin(found.id)].ipg)]
            missing = missing - set(storedIPGs.id)
            ipgs = ipgs[~ipgs.ipg.isin(storedIPGs.ipg)]
            more = self._fetch_from_sql(ipgs.id.unique().tolist())
            moreids = more.melt(id_vars=["assembly"], value_vars=self.column, var_name="type", value_name="id")
            stored = pd.concat([stored,more])
            found = pd.concat([found,moreids], ignore_index=True)
            missing = missing - set(found.id)
            notinipgs = missing - set(ipgs.id)
            if len(notinipgs) > 0:
                self._add_to_missing(notinipgs,np.NaN,"No IPG")
                missing = missing - notinipgs
            drop = set()
            for lost in missing:
                a = ipgs.loc[ipgs.id == lost,'assembly']
                if len(a):
                    a = a.iloc[0]
                    self._add_to_missing(lost, a, "Not in local storage")
                drop.add(lost)
            missing = missing - drop

        # Verify all possible IDs
        if len(missing) > 0:
            self._add_to_missing(missing,np.NaN,"Not found in storage")
        if len(stored) == 0:
            return seqrecords_to_dataframe([])
        return stored

    def fetchone(self, proteins, ipgs=None):
        """
        Fetch each gene neighborhood iteratively.

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db.ncbi import entrez
          >>> from rotifer.db.ncbi import ftp
          >>> ic = ncbi.IPGCursor(batch_size=1)
          >>> gnc = ftp.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.fetchall(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        Generator of rotifer.genome.data.NeighborhoodDF
        """
        found = self.__getitem__(proteins, ipgs=ipgs)
        for bid, block in found.groupby('block_id'):
            yield block.copy()

    def fetchall(self, proteins, ipgs=None):
        """
        Fetch all gene neighborhoods in a single dataframe.

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to identify identical
          sequences using NCBI's IPG reports. Example:

          >>> from rotifer.db.ncbi import entrez
          >>> from rotifer.db.sql import sqlite3
          >>> ic = ncbi.IPGCursor(batch_size=1)
          >>> gnc = sqlite3.GeneNeighborhoodCursor("mydb.db")
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.fetchall(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        return self.__getitem__(proteins, ipgs=ipgs)

    def stored(self, block):
        """
        Identify gene neighborhoods present in the SQLite3 storage.

        The block_id column is used to identify gene neighborhoods
        and exact matches in the SQL database are reported as True.

        Parameters
        ----------
        block: rotifer.genome.data.NeighborhoodDF
          Gene neighborhoods to check for presence in the database.

        Returns
        -------
        Pandas Series

        The Series elements are booleans: True if found in the
        storage, False otherwise.
        """
        sql = []
        if self._has_table("neighborhoods"):
            sql = pd.read_sql("SELECT DISTINCT block_id from neighborhoods", self._dbconn).block_id
        return block.block_id.isin(sql)

    def insert(self, block):
        """
        Add one or more new gene neighborhoods to the local storage.

        Parameters
        ----------
        block: rotifer.genome.data.NeighborhoodDF
          Gene neighborhood dataframe
        """
        block = block[~self.stored(block)]
        if len(block) > 0:
            block.to_sql("neighborhoods", self._dbconn, if_exists = 'append', index=False)
