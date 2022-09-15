import os
import sys
import types
import socket
import typing
import logging
import pandas as pd
from tqdm import tqdm
from ftplib import FTP
import rotifer
import rotifer.db.core
from rotifer import GlobalConfig
from rotifer.db.ncbi import NcbiConfig
from rotifer.db.ncbi import utils as rdnu
from rotifer.genome.utils import seqrecords_to_dataframe
logger = rotifer.logging.getLogger(__name__)

# Classes

class cursor():
    """
    This class represents a live connection to the NCBI FTP server.
    """
    def __init__(self, url=NcbiConfig['ftpserver'], tries=3, timeout=50, cache=GlobalConfig['cache']):
        """
        Parameters
        ----------
        url: string
          NCBI FTP server URL
        tries: integer
          Maximum number of attempts to download a file
        timeout: integer
          Maximum amount of time, in seconds, to wait 
          for a connection to the server
        cache: directory path
          Folder to store temporary files
        """
        self._error = 0
        self.url = url
        self.tries = tries
        self.timeout = timeout
        self.cache = cache
        self.connect()

    def connect(self):
        """
        Connect or reconnect to server.
        """
        import time
        self._error = 0
        attempt = 0
        while attempt < self.tries:
            try:
                if hasattr(self,"connection"):
                    self.connection.sendcmd("NOOP")
                else:
                    self.connection = FTP(self.url, timeout=self.timeout)
                    self.connection.login()
                break
            except:
                try:
                    self.connection = FTP(self.url, timeout=self.timeout)
                    self.connection.login()
                    self._error = 0
                    break
                except socket.timeout:
                    time.sleep(1)
                except TimeoutError:
                    time.sleep(1)
                self._error = 1
            attempt += 1

    # Load NCBI assembly reports
    def ftp_get(self, target, avoid_collision=False, outdir=None):
        '''
        Download a file from the NCBI's FTP site.

        Usage:
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.cursor(tries=3)
          localpath = ftp.ftp_get('genomes/README.txt')

        Parameters:
          target : string
            URL or relative path of the file
          avoid_collision : boolean, default False
            Avoid name collision by adding random suffixes
          outdir : string, default is the object's cache
            Output directory

        Returns:
          Path to the downloaded file
        '''
        from tempfile import NamedTemporaryFile
        from rotifer.db.ncbi import NcbiConfig
        logger.debug(f'downloading {NcbiConfig["ftpserver"]}{target}')

        # Create output directory, if necessary
        if not outdir:
            outdir = self.cache
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except:
                logger.error(f'failed to create download directory {outdir}')
                self._error = 1
                return None

        # Retrieve contents for each folder
        # Prepare local file handle
        if avoid_collision:
            parts = os.path.splitext(os.path.basename(target))
            prefix = parts[0] + '.' 
            suffix = None if parts[1] == '' else parts[1]
            outfh = NamedTemporaryFile(mode='w+b', suffix=suffix, prefix=prefix, dir=outdir, delete=False)
            outfile = outfh.name
        else:
            outfile = os.path.join(outdir, os.path.basename(target))
            outfh   = open(outfile,'wb')

        # To avoid problems with very long names, I change to the target
        # directory and, later, back to /
        p = target.replace('ftp://',"")
        p = p.replace(NcbiConfig['ftpserver'],"")
        p = p[1:] if p[0] == "/" else p
        p = p.split("/")
        self.connect()
        for i in list(range(len(p)-1)):
            self.connection.cwd(p[i])
        try:
            self.connection.retrbinary("RETR " + p[-1], outfh.write)
        except:
            logger.error(f'Unable to download {target}')
            self._error = 1
            return None
        self.connection.cwd("/")
        outfh.close()

        # Return pandas object
        logger.debug(f'download complete {NcbiConfig["ftpserver"]}{target}')
        return outfile

    # List files in ftp directory
    def ftp_ls(self, targets):
        '''
        List contents of an FTP directory.

        Usage:
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.cursor()
          contents = ftp.ftp_ls('genomes')

        Parameters:
          target : (list of) string(s)
            Path of one or more diretory(ies)

        Returns:
          Pandas DataFrame
        '''
        import pandas as pd

        # Process targets
        d = []
        if not (isinstance(targets,list) or isinstance(targets,tuple)):
            targets = [targets]
        self.connect()
        for target in targets:
            try:
                for x in self.connection.mlsd(target):
                    if x[0] == "." or x[0] == "..":
                        continue
                    x[1]["target"] = target
                    x[1]["name"] = x[0]
                    d.append(x[1])
            except:
                logger.error(f'''Could not retrive list for directory {target} at the NCBI's FTP site.''')
                continue
        if not d:
            self._error = 1
        d = pd.DataFrame(d)
        return d

    # Mimick opening of local files
    def ftp_open(self, target,  mode='rt', avoid_collision=True, delete=True):
        '''
        Open a file stored at the NCBI FTP site.

        Note:
          Compressed files will be uncompressed on the fly.

        Usage:
          # Open a genome's GBF file
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.cursor()
          path = "genomes/all/GCA/900/547/725/GCA_900547725.1_UMGS1014/"
          path += "GCA_900547725.1_UMGS1014_genomic.gbff.gz"
          fh = ftp.ftp_open(path, mode="rt")

        Parameters:
          target : string
            Relative path of the target file
          mode : string
            Set read mode ('r', 'rt', 'rb', etc) to open the file
          avoid_collision : boolean, boolean default True
            Avoid name collision with random suffixes
          delete : boolean, default True
            Whether files should be automatically removed when closed
            If set to False, the file will remain in the cache directory

        Returns:
          Open data stream (file handle-like) object
        '''
        from tempfile import _TemporaryFileWrapper
        import rotifer.core.functions as rcf
        self.connect()
        outfile = self.ftp_get(target, avoid_collision=avoid_collision)
        if self._error:
            return None
        else:
            return _TemporaryFileWrapper(rcf.open_compressed(outfile, mode), outfile, delete)

    def open_genome(self, accession, assembly_reports=None):
        """
        Open the GBFF file of a genome from the NCBI FTP site.

        Usage:
          # Just open
          
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.cursor()
          fh = ftp.open_genome("GCA_900547725.1")

          # Open and parse to a list of Bio.SeqRecord

          from rotifer.db.ncbi import ftp as ncbiftp
          from Bio import SeqIO
          ftp = ncbiftp.cursor()
          fh = ftp.open_genome("GCA_900547725.1")
          s = [ x for x in SeqIO.parse(fh, "genbank") ]

        Parameters:
          accession : string
            Genome accession
          assembly_reports: pandas DataFrame
            (Optional) NCBI's ASSEMBLY_REPORTS tables
            This dataframe can be loade using:

              from rotifer.db.ncbi as ncbi
              ar = ncbi.assembly_reports()

        Returns:
          Open data stream (file handle-like) object
        """
        import rotifer.core.functions as rcf

        # find genome and download
        path = self.genome_path(accession, assembly_reports=assembly_reports)
        if len(path) == 0:
           return None

        # Download checksum
        md5url = "/".join([path[0],"md5checksums.txt"])
        self.connect()
        for attempt in range(0,self.tries):
            md5 = self.ftp_open(md5url, mode='rt', avoid_collision=True, delete=True)
            if self._error:
                logger.error(f'''Could not fetch checksum for {accession}.''')
            else:
                try:
                    md5 = pd.read_csv(md5, sep=' +', names=['md5','filename'], engine="python")
                    md5 = md5[md5.filename.fillna("_").str.contains('_genomic.gbff.gz')]
                    md5 = md5.md5.iloc[0]
                    if md5:
                        break
                except:
                    logger.error(f'''Parsing of checksum for {accession} failed.''')
        if not md5:
            self._error = 1
            return None

        # Download genome
        gz = None
        self.connect()
        for attempt in range(0,self.tries):
            gz = self.ftp_open("/".join(path),  mode='rt', avoid_collision=True, delete=True)
            if self._error:
                logger.error(f'''Could not open GBFF for {accession}.''')
            else:
                md5gz = rcf.md5(gz.name)
                if md5 == md5gz:
                    break

        # Return file object
        return gz

    def genome_path(self, accession, assembly_reports=None):
        """
        Fetch the path of a genome at the NCBI FTP site.

        Usage:
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.cursor()
          path = ftp.genome_path("GCA_900547725.1")
          print("/".join(path))

        Parameters:
          accession : string
            Genome accession
          assembly_reports: pandas DataFrame
            (Optional) NCBI's ASSEMBLY_REPORTS tables
            This dataframe can be loade using:

              from rotifer.db.ncbi as ncbi
              ar = ncbi.assembly_reports()
        
        Returns:
          A tuple of two strings, empty when the genome is not found
        """
        from rotifer.db.ncbi import NcbiConfig
        path = ()

        # Extract genome path from assembly reports
        if isinstance(assembly_reports, pd.DataFrame) and not assembly_reports.empty:
            path = assembly_reports.query(f'assembly == "{accession}"')
            if not path.empty:
                path = path.ftp_path.iloc[0]
                path = path.replace(f'ftp://{NcbiConfig["ftpserver"]}','')
                path = (path,os.path.basename(path) + "_genomic.gbff.gz")

        # Retrieve genome path for newest version
        if len(path) == 0:
            path = accession[0:accession.find(".")].replace("_","")
            path = "/".join([ path[i : i + 3] for i in range(0, len(path), 3) ])
            path = f'/genomes/all/{path}'
            path = self.ftp_ls(path)
            if path.empty:
                return ()
            path = path.query(f'name.str.contains("{accession}")')
            path = path.sort_values(['name'], ascending=False).iloc[0]
            path = path.target + "/" + path['name']

            # Retrieve GBFF path
            if not path:
                return ()
            path = self.ftp_ls(path)
            if path.empty:
                return ()
            path = path[path['name'].str.contains(".gbff.gz")]
            if path.empty:
                return ()
            path = (path.target.iloc[0], path['name'].iloc[0])

        return path

    def genome_report(self, accession):
        """
        Fetch genome assembly reports.

        Usage:
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.cursor()
          contigs, assembly = ftp.genome_report("GCA_900547725.1")

        Parameters:
          accession : string
            Genome accession
        
        Returns:
          A tuple containing a Pandas DataFrame and a Pandas Series
          
          The Pandas DataFrame lists the assembly's contigs

          The Series contains the assembly properties and is 
          similar to the table from rotifer.db.ncbi.assemblies()

        """
        # Column names
        arcolumn = f"""                Assembly name : assembly_name
                                       Organism name : organism_name
                                             Isolate : isolate
                                               Taxid : taxid
                                           BioSample : biosample
                                          BioProject : bioproject
                                           Submitter : submitter
                                                Date : submission_date
                                       Assembly type : assembly_type
                                        Release type : release_type
                                      Assembly level : assembly_level
                               Genome representation : representative
                                         WGS project : wgs
                                     Assembly method : assembly_method
                                     Genome coverage : genome_coverage
                               Sequencing technology : sequencing
                                     RefSeq category : refseq_category
                          GenBank assembly accession : genbank
                           RefSeq assembly accession : refseq
                                Excluded from RefSeq : excluded_from_refseq
    RefSeq assembly and GenBank assemblies identical : identical""".split("\n")
        arcolumn = [ x.strip().split(" : ") for x in arcolumn ]
        arcolumn = { x[0]:x[1] for x in arcolumn }

        # Opening data file
        ar = [['column','value'], ['assembly', accession]]
        sc = []
        path = self.genome_path(accession)
        report = self.ftp_ls(path[0])
        report = report[report.name.str.contains("_assembly_report.txt")]
        report = path[0] + "/" + report.name.iloc[0]
        report = self.ftp_open(report)

        inar = True
        for row in report:
            row = row.strip()
            if row == "#" or row[0:2] == "##":
                inar = False
            elif row[0:15] == "# Sequence-Name":
                sc.append(['assembly'] + row[2:].split("\t"))
            elif inar and row[0:2] == "# ":
                ar.append(row[2:].split(":"))
            elif row[0] != '#':
                sc.append([accession] + row.split("\t"))

        ar = pd.DataFrame(ar[1:], columns=ar[0])
        ar.value = ar.value.str.strip()
        ar.column = ar.column.replace(arcolumn)
        sc = pd.DataFrame(sc[1:], columns=sc[0])
        ar = ar.set_index("column")

        return sc, ar

#class GenomeCursor(BaseCursor):
class GenomeCursor(rotifer.db.core.SimpleParallelProcessCursor):
    """
    Fetch genome sequences from the NCBI FTP site.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
    >>> import rotifer.db.ncbi as ncbi
    >>> ar = ncbi.assemblies(taxonomy=True)
    >>> ar = ar.query('superkingdom != "Eukaryota"')
    >>> g = ar.assembly.sample(10)
    >>> from rotifer.db.ncbi import ftp
    >>> gc = ftp.GenomeCursor(g, assembly_reports=ar)
    >>> genomes = gc.fetch_all()

    Parameters
    ----------
    assembly_reports: Pandas Dataframe
      Table of genomes
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
            assembly_reports=None,
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
            cache=GlobalConfig['cache'],
        ):
        super().__init__(
            progress=progress,
            tries=tries,
            batch_size=batch_size,
            threads=threads
        )
        self.assembly_reports = assembly_reports
        self.cache = cache

    def _getids(self,obj):
        return {obj.assembly}

    def _fetcher(self, accession, assembly_reports=None):
        import rotifer.db.ncbi.ftp as ncbiftp
        gfc = ncbiftp.cursor(tries=1, cache=self.cache)
        return gfc.open_genome(accession, assembly_reports=self.assembly_reports)

    def _parser(self, stream, accession):
        from Bio import SeqIO
        stack = []
        for s in SeqIO.parse(stream,"genbank"):
            s.assembly = accession
            stack.append(s)
        return stack

    def __getitem__(self, accession):
        """
        Download, parse and load a genome.

        Returns
        -------
        List of Bio.SeqRecord objects
        """
        objlist = []
        for attempt in range(0,self.tries):
            # Download and/or open data file
            try:
                stream = self._fetcher(accession, assembly_reports=self.assembly_reports)
            except RuntimeError:
                logger.error(f'Runtime error: '+str(sys.exc_info()[1]))
                continue
            except:
                logger.debug(f"Failed to download genome {accession}: {sys.exc_info()}")
                continue

            # Use parser and process results
            try:
                objlist = self._parser(stream, accession)
                break
            except:
                logger.debug(f"Failed to parse genome {accession}: {sys.exc_info()}")

        if len(objlist) == 0:
            logger.warn(f'Could not download genome {accession}')
            self.missing.add(accession)

        # Return data
        return objlist

    def _splitter(self, accessions):
        from rotifer.devel.alpha.gian_func import chunks
        size = self.batch_size
        if size == None or size == 0:
            size = max(int(len(accessions) / self.threads),1)
        return chunks(list(accessions), size)

    def _worker(self,accessions):
        stack = []
        for accession in accessions:
            objlist = self[accession]
            if len(objlist) != 0:
                stack.extend(objlist)
        return stack

class GenomeFeaturesCursor(GenomeCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
    >>> import rotifer.db.ncbi as ncbi
    >>> ar = ncbi.assemblies(taxonomy=True)
    >>> ar = ar.query('superkingdom != "Eukaryota"')
    >>> g = ar.assembly.sample(10)
    >>> from rotifer.db.ncbi import ftp
    >>> gfc = ftp.GenomeFeaturesCursor(g, assembly_reports=ar)
    >>> df = gfc.fetch_all()

    Parameters
    ----------
    assembly_reports: Pandas Dataframe
      Table of genomes
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
            assembly_reports=None,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
            cache=GlobalConfig['cache'],
        ):
        super().__init__(
            assembly_reports=assembly_reports,
            progress=progress,
            tries=tries,
            batch_size=batch_size,
            threads=threads,
            cache=cache,
        )
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

    def _getids(self,obj):
        return set(obj.assembly)

    def __getitem__(self, accession):
        item = super().__getitem__(accession)
        if len(item) == 0:
            item = seqrecords_to_dataframe(item)
        return item

    def _parser(self, stream, accession):
        from Bio import SeqIO
        stream = SeqIO.parse(stream,"genbank")
        stream = seqrecords_to_dataframe(
            stream,
            exclude_type = self.exclude_type,
            autopid = self.autopid,
            assembly = accession,
            codontable = self.codontable,
        )
        return stream

    def _worker(self, accessions):
        stack = []
        for accession in accessions:
            df = self.__getitem__(accessions)
            if len(df) != 0:
                stack.append(df)
        return stack

    def fetch_all(self, accessions):
        """
        Fetch genomes.

        Parameters
        ----------
        accessions: list of strings
          NCBI genomes accessions

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        stack = []
        for df in self.fetch_each(accessions):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

class GeneNeighborhoodCursor(GenomeFeaturesCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
    >>> from rotifer.db.ncbi import ftp
    >>> gfc = ftp.GeneNeighborhoodCursor(progress=True)
    >>> df = gfc.fetch_all(["EEE9598493.1"])

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
    min_block_id : int
      Starting number for block_ids, useful if calling this method
      multiple times
    assembly_reports: Pandas Dataframe
     Table of genomes
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
            column = 'pid',
            before = 7,
            after = 7,
            min_block_distance = 0,
            strand = None,
            fttype = 'same',
            min_block_id = 1,
            assembly_reports=None,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
            cache=GlobalConfig['cache'],
        ):
        super().__init__(
            assembly_reports = assembly_reports,
            exclude_type = exclude_type,
            autopid = autopid,
            codontable = codontable,
            progress = progress,
            tries = tries,
            batch_size = batch_size,
            threads = threads,
            cache = cache,
        )
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
        self.min_block_id = min_block_id

    def __getitem__(self, protein, ipgs=None):
        """
        Find gene neighborhoods in a genome.

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
            ipgs = ic.fetch_all(protein)
        ipgs = ipgs[ipgs.id.isin(ipgs[ipgs.pid.isin(protein) | ipgs.representative.isin(protein)].id)]
        best = rdnu.best_ipgs(ipgs)
        best = best[best.assembly.notna()]
        ipgs = ipgs[ipgs.assembly.isin(best.assembly)]
        if len(ipgs) == 0:
            self.missing.union(protein)
            return objlist

        # Identify DNA data
        assemblies, nucleotides = rdnu.ipgs_to_dicts(ipgs, best=False, full=False)

        # Download and parse
        objlist = []
        for accession in assemblies.keys():
            obj = None
            for attempt in range(0,self.tries):
                # Download and open data file
                try:
                    stream = self._fetcher(accession, assembly_reports=self.assembly_reports)
                except RuntimeError:
                    logger.error(f'Runtime error: '+str(sys.exc_info()[1]))
                    continue
                except:
                    logger.debug(f"Failed to download genome {accession}: {sys.exc_info()}")
                    continue

                # Use parser to process results
                try:
                    obj = self._parser(stream, accession, assemblies[accession])
                    break
                except:
                    logger.debug(f"Failed to parse genome {accession}: {sys.exc_info()}")

            if isinstance(obj, types.NoneType):
                logger.warn(f'Could not download genome {accession}')
            elif len(obj) == 0:
                logger.warn(f'No anchors in genome {accession}')
            else:
                objlist.append(obj)

        # No data?
        if len(objlist) == 0:
            return seqrecords_to_dataframe([])

        # Concatenate and evaluate
        objlist = pd.concat(objlist, ignore_index=True)
        found = set(objlist.pid).union(set(objlist.replaced))
        found = set(protein).intersection(found)
        self.missing = self.missing.union(set(protein) - found)

        # Return data
        return objlist

    def _parser(self, stream, accession, proteins):
        stream = super()._parser(stream, accession)
        stream = stream.neighbors(
            stream[self.column].isin(proteins.keys()),
            before = self.before,
            after = self.after,
            min_block_distance = self.min_block_distance,
            strand = self.strand,
            fttype = self.fttype,
            min_block_id = self.min_block_id
        )
        stream['replaced'] = stream.pid.replace(proteins)
        return stream

    def _worker(self, chunk):
        df = self.__getitem__(*chunk)
        return [ x[1] for x in df.groupby('block_id') ]

    def _splitter(self, ipgs):
        for x, y in ipgs.groupby('assembly'):
            proteins = set(y.pid).union(y.representative)
            yield (proteins, y.copy())

    def fetch_each(self, proteins, ipgs=None):
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
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries, batch_size=self.batch_size, threads=self.threads)
            ipgs = ic.fetch_all(list(proteins))
        ipgs = ipgs[ipgs.pid.isin(proteins) | ipgs.representative.isin(proteins)]
        if len(ipgs) == 0:
            self.missing = self.missing.union(proteins)
            yield seqrecords_to_dataframe([])
        assemblies = rdnu.best_ipgs(ipgs)
        assemblies = assemblies[assemblies.assembly.notna()]
        assemblies = ipgs[ipgs.assembly.isin(assemblies.assembly)]

        # Split jobs and execute
        last_block_id = 1
        todo = set(assemblies.assembly.unique())
        self.missing = self.missing.union(proteins)
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                logger.info(f'Downloading {len(todo)} genomes...')
                p = tqdm(total=len(todo), initial=0)
            tasks = []
            for chunk in self._splitter(assemblies):
                tasks.append(executor.submit(self._worker, chunk))
            for x in as_completed(tasks):
                for obj in x.result():
                    found = self._getids(obj)
                    done = todo.intersection(found)
                    ok = set(obj.pid)
                    if 'replaced' in obj:
                        ok = ok.union(obj.replaced)
                    self.missing = self.missing - ok
                    if len(obj) > 0:
                        obj.block_id = last_block_id
                        last_block_id += 1
                    if self.progress and len(done) > 0:
                        p.update(len(done))
                    todo = todo - done
                    yield obj

    def fetch_all(self, proteins, ipgs=None):
        """
        Fetch genomes.

        Parameters
        ----------
        proteins: list of strings
          NCBI protein identifiers

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        stack = []
        for df in self.fetch_each(proteins, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

# Is this library being used as a script?
if __name__ == '__main__':
    pass
