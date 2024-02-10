__doc__ = """
Rotifer connections to SQL databases
====================================
"""

# Dependencies
import re
import os
import sys
import uuid
import types
import typing
import sqlite3
import numpy as np
import pandas as pd
from tqdm import tqdm

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.ncbi.utils as rdnu
import rotifer.devel.beta.sequence as rdbs
from rotifer.core.functions import loadConfig
from rotifer.genome.data import NeighborhoodDF
logger = rotifer.logging.getLogger(__name__)
config = loadConfig(__name__, defaults = {})

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, rotifer.db.core.BaseCursor):
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
    >>> gnc = rdss.GeneNeighborhoodCursor("genomes.sqlite3")
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
      Path to a SQLite3 database file
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
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):

        super().__init__(progress = progress, *args, **kwargs)
        self.tries = tries
        self.batch_size = batch_size
        self.threads = threads
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
        self.path = path
        self.replace = replace
        self.cache = cache
        if os.path.exists(self.path) and self.replace:
            os.remove(self.path)
        self._dbconn = sqlite3.connect(self.path)
        self.uuid = str(uuid.uuid4())

    def __getitem__(self, accession, ipgs=None):
        """
        Dictionary-like access to gene neighbors.
        """
        if not self.has_table('features'):
            return NeighborhoodDF()
        if not isinstance(accession, typing.Iterable) or isinstance(accession,str):
            accession = [accession]

        # Register queries in the database and search neighborhoods
        self.submit(accession)
        sqlquery = f"""
            SELECT f3.nucleotide, f3.start, f3.end, f3.strand, t.block_id,
                   CASE WHEN t.ids LIKE "%" || f3.pid || "%" THEN 1 ELSE 0 END as query,
                   f3.pid, f3.type, f3.plen, f3.locus, f3.seq_type, f3.assembly, gene, f3.origin,
                   f3.topology, f3.product, f3.organism, f3.lineage, f3.classification,
                   f3.feature_order, f3.internal_id, f3.pid as replaced
            FROM (
                SELECT assembly, nucleotide, type, block_id, min(idup) as idup, max(iddown) as iddown, group_concat(id,char(1)) as ids
                FROM (
                    SELECT *, SUM(nooverlap) OVER (ORDER BY assembly, nucleotide, idup, iddown) as block_id
                    FROM (
                        SELECT *,
                               CASE WHEN 
                                 nucleotide = LAG(nucleotide) OVER (ORDER BY assembly, nucleotide, idup, iddown)
                                 and idup - LAG(iddown) OVER (ORDER BY assembly, nucleotide, idup, iddown) <= {self.min_block_distance}
                                THEN 0
                                ELSE 1
                               END AS nooverlap
                        FROM (
                            SELECT q.id, f1.assembly, f1.nucleotide, f1.type,
                                   f1.feature_order - {self.before} as foup, f1.feature_order + {self.after} as fodown,
                                   min(f2.internal_id) as idup, max(f2.internal_id) as iddown
                            FROM queries as q
                             inner join features as f1 on (q.id = f1.{self.column})
                             inner join features as f2 on (
                                f1.assembly = f2.assembly
                                and f1.nucleotide = f2.nucleotide
                                and f1.type == f2.type
                                and f2.feature_order >= foup
                                and f2.feature_order <= fodown
                             )
                            WHERE q.uuid = '{self.uuid}'
                            GROUP BY q.id, f1.assembly, f1.nucleotide, f1.type, foup, fodown
                            ORDER BY f1.assembly, f1.nucleotide, idup, iddown
                        ) as v
                    ) as w
                ) as z
                GROUP BY assembly, nucleotide, type, block_id
            ) as t
            inner join features as f3 on (
               t.assembly = f3.assembly 
               and t.nucleotide = f3.nucleotide
               and f3.internal_id >= idup 
               and f3.internal_id <= iddown
            )
            WHERE f3.type NOT IN ('{"','".join(self.exclude_type)}')
            ORDER BY f3.assembly, f3.nucleotide, f3.block_id, f3.start, f3.end
        """
        df = pd.read_sql(sqlquery, self._dbconn)
        self.cleanup()

        # Restrict results to nucleotides found in the IPG reports
        if not (isinstance(ipgs,type(None)) or ipgs.empty):
            df = df[df.nucleotide.isin(ipgs.nucleotide)]

        # Find missing entries, if any
        missing = set(accession).difference(self.getids(df, ipgs=ipgs))
        if len(missing):
            self.update_missing(missing, error=f'Entry not found in SQLite3 database {self.path}', retry=True)

        return NeighborhoodDF(df)

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

          >>> from rotifer.db.sql import sqlite3 as dns
          >>> gnc = rdss.GeneNeighborhoodCursor(progress=True)
          >>> for n in gnc.fetchone(['WP_063732599.1']):
                print(n.groupby('nucleotide').block_id.nunique())

        Returns
        -------
        Generator of rotifer.genome.data.NeighborhoodDF
        """
        if not isinstance(proteins,typing.Iterable) or isinstance(proteins,str):
            proteins = [proteins]
        if self.progress:
            logger.warn(f'Searching {len(proteins)} protein(s) in SQLite3 database at {self.path}')
            p = tqdm(total=len(proteins), initial=0)
        found = self.__getitem__(proteins, ipgs=ipgs)
        for bid, block in found.groupby('block_id'):
            done = self.getids(block, ipgs)
            done = proteins.intersection(done)
            if self.progress and len(done) > 0:
                p.update(len(done))
            yield block.copy()

    def fetchall(self, ids, ipgs=None):
        """
        Fetch all gene neighborhoods at once.

        Parameters
        ----------
        proteins: list of strings
          Database identifiers.
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI several times. Example:

          >>> from rotifer.db.sql import sqlite3 as rdss
          >>> gnc = rdss.GeneNeighborhoodCursor()
          >>> n = gnc.fetchall(['WP_063732599.1'])

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        return self.__getitem__(ids, ipgs=ipgs)

    def insert(self, data):
        """
        Store genome annotation data in the SQLite3 database.

        Parameters
        ----------
        data: rotifer.genome.data.NeighborhoodDF
          Gene neighborhood dataframe
        """
        data = data[~self.stored(data)]
        if len(data) > 0:
            data.to_sql("features", self._dbconn, if_exists = 'append', index=False)

    def stored(self, data, column='block_id', table='features'):
        """
        Verify if data is stored in the SQLite3 database.

        Parameters
        ----------
        data: string, list, pandas series, dataframe or derived
              classes such as rotifer.genome.data.NeighborhoodDF
        column: (list of) string
          Column(s) to search.
        table: string
          Name of the table to search for matches.

        Returns
        -------
        Pandas Series

        The Series elements are booleans, i.e. True if found in the
        storage, False otherwise. 
        """
        ret = pd.Series([ False for x in range(1,len(data)) ])
        if not self.has_table(table):
            return ret
        if isinstance(data, str):
            data = [data]
        if isinstance(data, list):
            data = pd.Series(data, name='input')
        for col in column:
            inStore = pd.read_sql(f"""SELECT DISTINCT {col} from {table}""", self._dbconn)[col]
            if isinstance(data, pd.DataFrame):
                ret = ret | data[col].isin(inStore)
            else:
                ret = ret | data.isin(inStore)
        return ret

    def has_table(self, name):
        """
        Find whether a table exists in the database.
        """
        sql = self._dbconn.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{name}'").fetchall()
        return len(sql) > 0

    def schema(self):
        """
        Show table definitions.

        Returns
        -------
        String
        """
        return self._dbconn.execute("""SELECT sql FROM sqlite_schema;""").fetchall()[0][0]

    def submit(self, accessions):
        """
        Send data to temporary tables for use in subsequent searches.
        
        Note
        ----
        Submitted data will be overwritten in the next submit() call.

        Parameters
        ----------
        accessions:
         Any sqlite3 supported values, such as strings, integers or float.
        """
        ids = set(accessions)
        cursor = self._dbconn.cursor()
        cursor.execute('CREATE TEMPORARY TABLE IF NOT EXISTS queries (id, uuid)')
        self.cleanup()
        cursor.executemany("INSERT INTO queries VALUES (?,?)",[ (x,self.uuid) for x in ids ])
        self._dbconn.commit()

    def cleanup(self):
        """
        Remove data from temporary tables.
        """
        self._dbconn.execute(f"DELETE FROM queries WHERE uuid = '{self.uuid}'")
        self._dbconn.commit()
