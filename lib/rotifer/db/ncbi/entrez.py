# Import external modules
import os
import sys
import types
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
from rotifer.db.ncbi import utils as rdnu
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
    >>> seqrec = eutils.fetchall("YP_009724395.1")

    Fetch several nucleotide entries
    >>> import sys
    >>> from Bio import SeqIO
    >>> import rotifer.db.ncbi as ncbi
    >>> eutils = ncbi.SequenceCursor(database="nucleotide")
    >>> query = ['CP084314.1', 'NC_019757.1', 'AAHROG010000026.1']
    >>> for seqrec in eutils.fetchone(query):
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

    def getids(self, obj):
        return {obj.id}

    def parser(self, stream, accession):
        stack = []
        for s in SeqIO.parse(stream, self._format):
            stack.append(s)
        return stack

    def fetcher(self, accession):
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
                    stream = self.fetcher(acc)
                except RuntimeError:
                    logger.error(f'Runtime error: '+str(sys.exc_info()[1]))
                    continue
                except:
                    logger.debug(f"Efetch failed for {accession}: {sys.exc_info()}")
                    continue
                try:
                    objlist = self.parser(stream, acc)
                except:
                    logger.debug(f"Parser failed for {accession}: {sys.exc_info()}")
                    continue
                for obj in objlist:
                    found = found.union(self.getids(obj))
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

    def worker(self,accessions):
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

    def fetchone(self,accessions):
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
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                p = tqdm(total=len(accessions), initial=0)
            tasks = []
            for chunk in chunks(list(accessions), size):
                tasks.append(executor.submit(self.worker, chunk))
            todo = accessions
            for x in as_completed(tasks):
                for obj in x.result():
                    found = self.getids(obj)
                    done = todo.intersection(found)
                    self.missing = self.missing - done
                    if self.progress and len(done) > 0:
                        p.update(len(done))
                    todo = todo - done
                    yield obj

    def fetchall(self, accessions):
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
        return list(self.fetchone(accessions))

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

    def getids(self,obj):
        return set(obj.pid).union(obj.representative)

    def parser(self, stream, accession):
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
            rep = ipg.query('is_query == 1').drop_duplicates('id', keep='first')
            rep = rep.set_index('id').pid.to_dict()
            ipg['representative'] = ipg['id'].map(rep)
            # Remove IPGs with no known query
            ipg = ipg[ipg.representative.notna()]
        else: # One query
            ipg['representative'] = accession

        # Register all accessions found and return sliced DataFrame
        return [ x[1].copy() for x in ipg.groupby('id') ]

    def fetchone(self,accessions):
        seen = set()
        for ipg in super().fetchone(accessions):
            ipg = ipg[~ipg.id.isin(seen)]
            if len(ipg) == 0:
                continue
            seen = seen.union(ipg.id)
            yield ipg

    def fetchall(self, accessions):
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame(columns = self._columns + self._added_columns)
        return df

class TaxonomyCursor(SequenceCursor):
    def __init__(self,progress=False, tries=3, batch_size=None, threads=10):
        super().__init__(database="taxonomy", progress=progress, tries=tries, batch_size=batch_size, threads=threads)
        self._rettype = "full"
        self._retmode = 'xml'

    def getids(self,obj):
        other = obj.alternative_taxids.str.split(",").explode().dropna()
        return set(obj.taxid).union(other)

    def parser(self, stream, accession):
        import Bio
        from Bio import Entrez
        taxdf = Entrez.parse(stream)
        taxdf = pd.DataFrame(Entrez.parse(stream))
        cols = [ taxdf[x].dropna().map(lambda x: isinstance(x,str)).all() for x in taxdf.columns.to_list() ]
        cols.append(True)
        if "AkaTaxIds" in taxdf.columns:
            taxdf["alternative_taxids"] = taxdf["AkaTaxIds"].fillna("").map(lambda x: ",".join(x)).replace("",np.nan)
        else:
            taxdf["alternative_taxids"] = np.nan
        taxdf = taxdf.loc[:,cols].applymap(lambda x: str(x))
        taxdf['superkingdom'] = taxdf.Lineage.str.replace("cellular organisms; ","").str.split("; ", expand=True)[0]
        taxdf.rename({'Lineage':'classification', 'TaxId':'taxid'}, axis=1, inplace=1)
        return [taxdf]

    def fetchall(self, accessions):
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = pd.DataFrame()
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

    def getids(self,obj):
        return set(obj.nucleotide)

    def parser(self, stream, accession):
        stream = SeqIO.parse(stream, self._format)
        stream = seqrecords_to_dataframe(stream, exclude_type=self.exclude_type, autopid=self.autopid, assembly=self.assembly, codontable=self.codontable)
        stream = [ x[1].copy() for x in stream.groupby('nucleotide') ]
        return stream

    def fetchall(self, accessions):
        df = list(self.fetchone(accessions))
        if len(df) > 0:
            df = pd.concat(df, ignore_index=True)
        else:
            df = seqrecords_to_dataframe([])

class GeneNeighborhoodCursor(NucleotideFeaturesCursor):
    """
    Fetch gene neighbors from nucleotide annotation.

    Usage
    -----
    >>> from rotifer.db.ncbi import ftp
    >>> gfc = ftp.GeneNeighborhoodCursor(progress=True)
    >>> df = gfc.fetchall(["EEE9598493.1"])

    Parameters
    ----------
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
      If set to True, neighborhood data for eukaryotic nucleotide sequences
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
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
        ):
        super().__init__(
            exclude_type = exclude_type,
            autopid = autopid,
            codontable = codontable,
            progress = progress,
            tries = tries,
            batch_size = batch_size,
            threads = threads,
        )
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
        self.eukaryotes = eukaryotes
        self.missing = pd.DataFrame(columns=["noipgs","eukaryote","assembly","error","class"])

    def _pids(self, obj):
        ids = obj.melt(id_vars=['nucleotide'], value_vars=['pid','replaced'], value_name='pid', var_name='type')
        ids = set(ids.dropna().pid)
        return ids

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

    def __getitem__(self, protein, ipgs=None):
        """
        Find gene neighborhoods in a nucleotide sequence.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        objlist = seqrecords_to_dataframe([])
        if not isinstance(protein,typing.Iterable) or isinstance(protein,str):
            protein = [protein]

        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            ic = entrez.IPGCursor(progress=False, tries=self.tries, batch_size=self.batch_size, threads=self.threads)
            ipgs = ic.fetchall(protein)
        ipgs = ipgs[ipgs.id.isin(ipgs[ipgs.pid.isin(protein) | ipgs.representative.isin(protein)].id)]
        best = rdnu.best_ipgs(ipgs)
        best = best[best.nucleotide.notna()]
        ipgs = ipgs[ipgs.nucleotide.isin(best.nucleotide)]
        missing = set(protein) - set(ipgs.pid) - set(ipgs.representative)
        if missing:
            self._add_to_missing(missing,np.NaN,"No IPGs")
        if len(ipgs) == 0:
            return objlist

        # Identify DNA data
        nucleotides = ipgs.filter(['nucleotide','pid','representative'])
        nucleotides = nucleotides.drop_duplicates(ignore_index=True)
        nucleotides = nucleotides.groupby('nucleotide').apply(lambda x: x.set_index('pid').representative.to_dict())
        nucleotides = nucleotides.to_dict()
        #assemblies, nucleotides = rdnu.ipgs_to_dicts(ipgs)

        # Download and parse
        objlist = []
        for accession in nucleotides.keys():
            expected = set([ y for x in nucleotides[accession].items() for y in x ])

            obj = None
            for attempt in range(0,self.tries):
                # Download and open data file
                error = None
                stream = None
                try:
                    stream = self.fetcher(accession)
                except RuntimeError:
                    error = f'Runtime error: {sys.exc_info()[1]}'
                    logger.debug(error)
                    continue
                except ValueError:
                    error = f'Value error: {sys.exc_info()[1]}'
                    logger.debug(error)
                    break
                except:
                    error = f'Failed to download nucleotide {accession}: {sys.exc_info()[1]}'
                    logger.debug(error)
                    continue

                if isinstance(stream, types.NoneType):
                    self._add_to_missing(expected, accession, error)
                    continue

                # Use parser to process results
                try:
                    obj = self.parser(stream, accession, nucleotides[accession])
                    break
                except ValueError:
                    error = f'Value error: {sys.exc_info()[1]}'
                    logger.debug(error)
                    break
                except:
                    error = f"Failed to parse nucleotide {accession}: {str(sys.exc_info()[1])}"
                    logger.debug(error)

            if isinstance(obj, types.NoneType):
                self._add_to_missing(expected, accession, error)
            elif len(obj) == 0:
                error = f'No anchors in nucleotide sequence {accession}'
                self._add_to_missing(expected, accession, error)
            else:
                objlist.extend(obj)

        # No data?
        if len(objlist) == 0:
            return seqrecords_to_dataframe([])

        # Concatenate and evaluate
        objlist = pd.concat(objlist, ignore_index=True)

        # Return data
        if len(objlist) > 0:
            self.missing.drop(self._pids(objlist), axis=0, inplace=True, errors='ignore')
        return objlist

    def parser(self, stream, accession, proteins):
        stack = []
        for df in super().parser(stream, accession):
            taxonomy = df.classification.fillna("").iloc[0].split(";")
            if (not self.eukaryotes) and "Eukaryota" in taxonomy:
                raise ValueError(f"Eukaryotic nucleotide sequence {accession} ignored.")
            df = df.neighbors(
                df[self.column].isin(proteins.keys()),
                before = self.before,
                after = self.after,
                min_block_distance = self.min_block_distance,
                strand = self.strand,
                fttype = self.fttype,
            )
            df['replaced'] = df.pid.replace(proteins)
            stack.append(df)
        return stack

    def worker(self, chunk):
        result = []
        for args in chunk:
            df = self.__getitem__(*args)
            if len(df) == 0:
                continue
            for x in df.groupby('block_id'):
                result.append(x[1])
        return {"result":result,"missing":self.missing}

    def splitter(self, ipgs):
        size = self.batch_size
        if size == None or size == 0:
            size = max(int(ipgs.nucleotide.nunique()/self.threads),1)
        batch = []
        for x, y in ipgs.groupby('nucleotide'):
            proteins = set(y.pid).union(y.representative)
            batch.append((proteins, y.copy()))
        batch = [ batch[x:x+size] for x in range(0,len(batch),size) ]
        return batch

    def fetchone(self, proteins, ipgs=None):
        """
        Asynchronously fetch gene neighborhoods from NCBI.

        Parameters
        ----------
        proteins: list of strings
          NCBI protein identifiers

        Returns
        -------
        A generator for rotifer.genome.data.NeighborhoodDF objects
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Make sure no identifiers are used twice
        if not isinstance(proteins,typing.Iterable) or isinstance(proteins,str):
            proteins = [proteins]
        proteins = set(proteins)

        # Make sure we have IPGs
        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            if self.progress:
                logger.info(f'Downloading IPGs for {len(proteins)} proteins...')
            size = self.batch_size
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries, threads=self.threads)
            ipgs = ic.fetchall(list(proteins))
            if len(ic.missing):
                self._add_to_missing(ic.missing.index.to_list(), np.nan, "No IPGs for nucleotides at NCBI")
        ipgs = ipgs[ipgs.pid.isin(proteins) | ipgs.representative.isin(proteins)]
        if len(ipgs) == 0:
            self._add_to_missing(proteins,np.NaN,"No IPGs to match a nucleotide sequence")
            return [seqrecords_to_dataframe([])]
        nucleotides = rdnu.best_ipgs(ipgs)
        nucleotides = nucleotides[nucleotides.nucleotide.notna()]
        nucleotides = ipgs[ipgs.nucleotide.isin(nucleotides.nucleotide)]
        if len(nucleotides) == 0:
            return [seqrecords_to_dataframe([])]

        # Split jobs and execute
        todo = set(nucleotides.nucleotide.unique())
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                targets = set(nucleotides.pid).union(nucleotides.representative)
                targets = len(targets.intersection(proteins))
                logger.info(f'Downloading {len(todo)} nucleotides for {targets} proteins...')
                p = tqdm(total=len(todo), initial=0)
            tasks = []
            for chunk in self.splitter(nucleotides):
                tasks.append(executor.submit(self.worker, chunk))
            completed = set()
            for x in as_completed(tasks):
                data = x.result()
                for s in data['missing'].iterrows():
                    self.missing.loc[s[0]] = s[1]
                for obj in data['result']:
                    found = self.getids(obj)
                    done = todo.intersection(found)
                    if self.progress and len(done) > 0:
                        p.update(len(done))
                    todo = todo - done
                    completed.update(set(obj.pid.dropna()).union(obj.replaced.dropna()))
                    found = completed.intersection(self.missing.index)
                    if found:
                        self.missing.drop(found, axis=0, inplace=True)
                    yield obj

    def fetchall(self, proteins, ipgs=None):
        """
        Fetch all neighbors from nucleotide sequences.

        Parameters
        ----------
        proteins: list of strings
          NCBI protein identifiers

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        stack = []
        for df in self.fetchone(proteins, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

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
