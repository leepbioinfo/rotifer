# Import external modules
import os
import sys
import socket
import typing
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

# Import submodules
import rotifer
from rotifer import GlobalConfig
from rotifer.db.ncbi import NcbiConfig
from rotifer.core.functions import findDataFiles
from rotifer.core.functions import loadConfig
from rotifer.genome.utils import seqrecords_to_dataframe

# Load NCBI subclasses
logger = rotifer.logging.getLogger(__name__)

# Classes

class SequenceCursor:
    """
    Fetch sequences from NCBI using accession numbers and EUtilities.

    This class downloads and parses sequences in Genbank format, i.e.
    the most richly annotated format from NCBI.

    Usage
    -----
    Fetch a protein sequence
    >>> from rotifer.db.ncbi import entrez
    >>> eutils = entrez.SequenceCursor(database="protein")
    >>> seqrec = eutils.fetch_all("YP_009724395.1")

    Fetch several nucleotide entries
    >>> import sys
    >>> from Bio import SeqIO
    >>> import rotifer.db.ncbi as ncbi
    >>> eutils = ncbi.SequenceCursor(database="nucleotide")
    >>> query = ['CP084314.1', 'NC_019757.1', 'AAHROG010000026.1']
    >>> for seqrec in eutils.fetch_each(query):
    >>>     print(SeqIO.write(seqrec, sys.stdout, "genbank")

    Parameters
    ----------
    database: string, default 'protein'
        Valid NCBI sequence database
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download data
    sleep_between_tries: int, default 1
      Number of seconds to wait between download attempts
    batch_size: int, default 1
      Number of accessions per batch
    threads: integer, default 3
      Number of simultaneous threads to run

    """
    def __init__(
            self,
            database="nucleotide",
            progress=False,
            tries=3,
            batch_size=None,
            threads=10
            ):
        self.missing = set()
        self.database = database
        self.progress = progress
        self.tries = tries
        self.batch_size = batch_size
        self.threads = threads
        if self.threads > 3:
            if NcbiConfig['api_key']:
                if self.threads > 10:
                    self.threads = 10
            else:
                self.threads = 3

        # Private attributes (may be overloaded by children)
        self._rettype = 'gbwithparts'
        self._format = 'genbank'
        self._retmode = 'text'

    def _getids(self, obj):
        return {obj.id}

    def _parser(self, stream, accession):
        stack = []
        for s in SeqIO.parse(stream, self._format):
            stack.append(s)
        return stack

    def _fetcher(self, accession):
        from Bio import Entrez
        Entrez.email = NcbiConfig["email"]
        Entrez.api_key = NcbiConfig["api_key"]
        return Entrez.efetch(db=self.database, rettype=self._rettype, retmode=self._retmode, id=accession, max_tries=self.tries, sleep_between_tries=1)

    def __getitem__(self, accession):
        """
        Run EFetch to download NCBI data.

        Parameters
        ----------
        accession: string
          Comma separated list of NCBI database entries.

        Returns
        -------
        List of Bio.SeqRecord objects
        """
        query = set(accession.split(","))
        batch = [accession]
        stack = []
        found = set()
        for attempt in range(0,2):
            for acc in batch:
                try:
                    stream = self._fetcher(acc)
                except RuntimeError:
                    logger.error(f'Runtime error: '+str(sys.exc_info()[1]))
                    continue
                except:
                    logger.debug(f"Efetch failed for {accession}: {sys.exc_info()}")
                    continue
                try:
                    objlist = self._parser(stream, acc)
                except:
                    logger.debug(f"Parser failed for Efetch of {accession}: {sys.exc_info()}")
                    continue
                for obj in objlist:
                    found = found.union(self._getids(obj))
                    stack.append(obj)
            if found == query:
                break
            else:
                batch = query - found

        missed = query - found
        if len(missed) > 0:
            logger.debug(f'''Unable to fetch data from {self.database} database for accessions {missed}''')
            self.missing = self.missing.union(missed)

        if len(stack) == 1 and len(query) == 1:
            stack = stack[0]

        return stack

    def _worker(self,accessions):
        """
        Fetch several accessions.

        Notes
        -----
        * This method is executed in the main thread

        Returns
        -------
        A list of Bio.SeqRecord objects
        """
        stack = []
        for chunk in [ accessions[x:x+200] for x in range(0,len(accessions),200) ]:
            it = self.__getitem__(",".join(chunk))
            if not isinstance(it, list):
                it = [it]
            for obj in it:
                stack.append(obj)
        return stack

    def fetch_each(self,accessions):
        """
        Asynchronously fetch sequences from NCBI.

        Notes
        -----
        * This method is executed in the main thread
        * Input order is not preserved.

        Parameters
        ----------
        accessions: list of strings
          List of NCBI accession numbers

        Returns
        -------
        A generator for Bio.SeqRecord objects
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from rotifer.devel.alpha.gian_func import chunks

        # Process arguments
        if not isinstance(accessions,typing.Iterable) or isinstance(accessions,str):
            accessions = [accessions]
        accessions = set(accessions)
        size = self.batch_size
        if size == None or size == 0:
            size = max(int(len(accessions) / self.threads),1)

        # Split jobs and execute
        self.missing = self.missing.union(accessions)
        todo = len(self.missing)
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                p = tqdm(total=len(accessions), initial=0)
            tasks = []
            for chunk in chunks(list(accessions), size):
                tasks.append(executor.submit(self._worker, chunk))
            for x in as_completed(tasks):
                for obj in x.result():
                    found = self._getids(obj)
                    self.missing = self.missing - found
                    done = todo - len(self.missing)
                    if self.progress and done:
                        p.update(done)
                    todo = len(self.missing)
                    yield obj

    def fetch_all(self, accessions):
        """
        Fetch all sequences from NCBI (blocking method).
        Note: input order is not preserved.

        Parameters
        ----------
        accessions: list of strings
          NCBI genomes accessions

        Returns
        -------
        List of Bio.SeqRecord objects
        """
        return list(self.fetch_each(accessions))

class FastaCursor(SequenceCursor):
    def __init__(self, database="protein", progress=False, tries=3, batch_size=None, threads=10):
        super().__init__(database=database, progress=progress, tries=tries, batch_size=batch_size, threads=threads)
        self._rettype = "fasta"
        self._format = 'fasta'

class IPGCursor(SequenceCursor):
    def __init__(self,progress=False, tries=3, batch_size=None, threads=10):
        super().__init__(database="protein", progress=progress, tries=tries, batch_size=batch_size, threads=threads)
        self._rettype = "ipg"
        self._columns = ['id','ipg_source','nucleotide','start','stop','strand','pid','description','ipg_organism','strain','assembly']
        self._added_columns = ['order','is_query','representative']

    def _getids(self,obj):
        return set(obj.pid).union(obj.representative)

    def _parser(self, stream, accession):
        query = accession.split(",")
        ipg = pd.read_csv(stream, sep='\t', names=self._columns, header=0).drop_duplicates()

        # Make sure all IPG ids are numbers
        numeric = pd.to_numeric(ipg.id, errors="coerce")
        errors = numeric.isna()
        if errors.any():
            logger.debug(f'Errors in IPG for accessions {accession}:\n{ipg[errors].to_string()}')
        ipg.id = numeric
        ipg = ipg[~errors]
        if ipg.empty:
            logger.debug(f'After removing errors, accessions {accession} were found to have no IPGs! Ignoring...')
            return []

        # Register query proteins
        ipg['is_query'] = ipg.pid.isin(query).astype(int).to_list()

        # Register original order of the table's rows
        o = pd.Series(range(1, len(ipg) + 1))
        c = pd.Series(np.where(ipg.id != ipg.id.shift(1), o.values, pd.NA)).ffill()
        ipg['order'] = (o - c).values

        # Annotate representatives
        if len(query) > 1: # Many queries
            #  Register first query protein as representative
            rep = ipg.loc[ipg['is_query'] == 1,['id','pid']].drop_duplicates('id', keep='first')
            rep = rep.set_index('id').pid.to_dict()
            ipg['representative'] = ipg['id'].map(rep)
            ipg = ipg[ipg.representative.notna()]
        else: # One query
            ipg['representative'] = accession

        # Register all accessions found and return sliced DataFrame
        return [ x[1].copy() for x in ipg.groupby('id') ]

    def fetch_each(self,accessions):
        seen = set()
        for ipg in super().fetch_each(accessions):
            ipg = ipg[~ipg.id.isin(seen)]
            if len(ipg) == 0:
                continue
            seen = seen.union(ipg.id)
            yield ipg

    def fetch_all(self, accessions):
        df = list(self.fetch_each(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame(columns = self._columns + self._added_columns)
        return df

class NucleotideFeaturesCursor(SequenceCursor):
    def __init__(
            self,
            exclude_type = ['source','gene','mRNA'],
            autopid = False,
            assembly = None,
            codontable= 'Bacterial',
            progress = False,
            tries = 3,
            batch_size = None,
            threads = 10
        ):
        super().__init__(database='nucleotide', progress=progress, tries=tries, batch_size=batch_size, threads=threads)
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.assembly = assembly
        self.codontable = codontable

    def _getids(self,obj):
        return set(obj.nucleotide)

    def _parser(self, stream, accession):
        stream = SeqIO.parse(stream, self._format)
        stream = seqrecords_to_dataframe(stream, exclude_type=self.exclude_type, autopid=self.autopid, assembly=self.assembly, codontable=self.codontable)
        stream = [ x[1].copy() for x in stream.groupby('nucleotide') ]
        return stream

    def fetch_all(self, accessions):
        df = list(self.fetch_each(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = seqrecords_to_dataframe([])

class GeneNeighborhoodCursor(NucleotideFeaturesCursor):
    from rotifer.db.ncbi.entrez import IPGCursor

    def __init__(
            self,
            column = 'pid',
            before = 3,
            after = 3,
            min_block_distance = 0,
            strand = None,
            eukaryotes = False,
            fttype = 'same',
            save = None,
            replace = True,
            exclude_type = ['source','gene','mRNA'],
            autopid = False,
            assembly = None,
            codontable= 'Bacterial',
            min_block_id = 1,
            progress = False,
            tries = 3,
            batch_size = None,
            threads = 10
        ):
        super().__init__(exclude_type=exclude_type,autopid=autopid,assembly=assembly,codontable=codontable,progress=progress,tries=tries,batch_size=batch_size,threads=threads)
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.eukaryotes = eukaryotes
        self.fttype = fftype,
        self.save = save
        self.replace = replace

    def _getids(self,obj):
        return set(obj.pid).union(obj.replaced)

    def _parser(self, stream, accession):
        stream = SeqIO.parse(stream, self._format)
        stream = seqrecords_to_dataframe(stream, exclude_type=self.exclude_type, autopid=self.autopid, assembly=self.assembly, codontable=self.codontable)
        stream = stream.neighbors()
        stream = [ x[1].copy() for x in stream.groupby('nucleotide') ]
        return stream

    def fetch_each(self, accessions, ipgs=None):
        # Process IPGs
        accessions = set(accessions)
        ipgs = self.ipgs
        if isinstance(ipgs,types.NoneType):
            logger.info(f'Downloading IPGs for {len(accessions)} accessions')
            ic = IPGCursor(progress=progress, tries=tries, batch_size=batch_size, threads=threads)
            ipgs = ic.fetch_all(accessions)
        ipgs = ipgs[ipgs.pid.isin(accessions) | ipgs.representative.isin(accessions)].copy()

        # Process each nucleotide sequence
        acc = set(ipgs.pid).union(ipgs.representative)
        for s in super().fetch_each(ipgs.nucleotide):
            s = s.neighbors(s[self.column].isin(acc), before=self.before, after=self.after, min_block_distance=self.min_block_distance)

    def fetch_all(self, accessions):
        df = list(self.fetch_each(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = seqrecords_to_dataframe([])
        return df

def elink(accessions, dbfrom="protein", dbto="taxonomy", linkname=None):
    """
    Find related database entries via NCBI's EUtilities.

    Usage:
      # download from NCBI's FTP site
      from rotifer.db.ncbi import entrez
      eutils = entrez.cursor()
      a = eutils.elink("YP_009724395.1")

    Returns:
      Pandas DataFrame

    Parameters:
      accessions : string or list of strings
        NCBI accessions to search links for
      to: string
        Name of the target database
      linkname: string
        Type of link between dbfrom and dbto
        If not set, {dbfrom}_{dbto} is used
    """
    from Bio import Entrez
    Entrez.email = NcbiConfig["email"]
    Entrez.api_key = NcbiConfig["api_key"]

    # Fix input
    if not isinstance(accessions,list):
        accessions = [accessions]
    if not linkname:
        linkname = dbfrom + "_" + dbto

    data = []
    for acc in accessions:
        try:
            raw = list(Entrez.read(Entrez.elink(dbfrom=dbfrom, linkname=linkname, id=acc)))
        except:
            logger.info(f'Entrez.elink failed for accession {acc}, dbfrom: {dbfrom}, dbto: {dbto}. Error: '+str(sys.exc_info()[0]))
            continue
        for d in raw:
            for x in d["LinkSetDb"]:
                for y in x["Link"]:
                    data.append([acc, d["IdList"][0], d["DbFrom"], x["LinkName"], dbto, y["Id"]])
    data = pd.DataFrame(data, columns=["qacc", "quid", "dbfrom", "linkname", "dbto", "tuid"])

    return data
