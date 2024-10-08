import os
import sys
import types
import socket
import typing
import numpy as np
import pandas as pd
from tqdm import tqdm
from ftplib import FTP
from copy import deepcopy

import rotifer
import rotifer.db.core
import rotifer.db.parallel
from rotifer.db.ncbi import NcbiConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Configuration
_defaults = {
    'batch_size': 4,
    "maxgetitem": 1,
    "threads": 15,
}
config = loadConfig(__name__, defaults = _defaults)

# Classes

class connection():
    """
    This class represents a live connection to the NCBI FTP server.
    """
    def __init__(self, url=NcbiConfig['ftpserver'], tries=3, timeout=50, cache=rotifer.config['cache']):
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
                    break
                except socket.timeout:
                    time.sleep(1)
                except TimeoutError:
                    time.sleep(1)
            attempt += 1

    # Download files
    def ftp_get(self, target, avoid_collision=False, outdir=None):
        '''
        Download a file from the NCBI's FTP site.

        Usage:
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.connection(tries=3)
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

        # Create output directory, if necessary
        if not outdir:
            outdir = self.cache
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except:
                raise IOError(f'failed to create download directory {outdir}')

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
        p = p.replace(self.url,"")
        p = p[1:] if p[0] == "/" else p
        p = p.split("/")
        self.connect()
        for i in list(range(len(p)-1)):
            self.connection.cwd(p[i])
        try:
            self.connection.retrbinary("RETR " + p[-1], outfh.write)
        except:
            raise IOError(f'Unable to write stream to {target}')
        self.connection.cwd("/")
        outfh.close()

        # Return pandas object
        logger.debug(f'Download complete {self.url}/{target}')
        return outfile

    # List files in ftp directory
    def ftp_ls(self, targets):
        '''
        List contents of an FTP directory.

        Usage:
          from rotifer.db.ncbi import ftp as ncbiftp
          ftp = ncbiftp.connection()
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
                raise FileNotFoundError(f'''Could not retrieve list for directory {target} at the NCBI's FTP site.''')
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
          ftp = ncbiftp.connection()
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
        return _TemporaryFileWrapper(rcf.open_compressed(outfile, mode), outfile, delete)

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Fetch genome sequences from the NCBI FTP site.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
    >>> from rotifer.db.ncbi import ftp
    >>> gc = ftp.GenomeCursor(g)
    >>> genomes = gc.fetchall()

    Parameters
    ----------
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
            progress=True,
            tries=3,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            timeout=10,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        self.timeout = timeout
        self.cache = cache

    def open_genome(self, accession, assembly_reports=None):
        """
        Open the GBFF file of a genome from the NCBI FTP site.

        Usage:
          # Just open
          
          from rotifer.db.ncbi import ftp
          gc = ftp.GenomeCursor()
          fh = gc.open_genome("GCA_900547725.1")

          # Open and parse to a list of Bio.SeqRecord

          from rotifer.db.ncbi import ftp
          from Bio import SeqIO
          gc = ftp.GenomeCursor()
          fh = gc.open_genome("GCA_900547725.1")
          s = [ x for x in SeqIO.parse(fh, "genbank") ]
          fh.close()

        Parameters:
          accession : string
            Genome accession
          assembly_reports: pandas DataFrame
            (Optional) NCBI's ASSEMBLY_REPORTS tables
            This dataframe can be loaded using:

              from rotifer.db.ncbi as ncbi
              ar = ncbi.assemblies()

        Returns:
          Open data stream (file handle-like) object
        """
        import rotifer.core.functions as rcf
        from rotifer.db.ncbi import ftp as ncbiftp
        ftp = ncbiftp.connection(tries=self.tries, timeout=self.timeout, cache=self.cache)

        # find genome and download
        path = self.genome_path(accession, assembly_reports=assembly_reports)
        if len(path) == 0:
           return None

        # Download checksum
        md5url = "/".join([path[0],"md5checksums.txt"])
        ftp.connect()
        for attempt in range(0,self.tries):
            md5 = ftp.ftp_open(md5url, mode='rt', avoid_collision=True, delete=True)
            try:
                md5 = pd.read_csv(md5, sep=' +', names=['md5','filename'], engine="python")
                md5 = md5[md5.filename.fillna("_").str.contains('_genomic.gbff.gz')]
                md5 = md5.md5.iloc[0]
                if md5:
                    break
            except:
                raise IOError(f'''Parsing of checksum for {accession} failed.''')
        if not md5:
            return None

        # Download genome
        gz = None
        ftp.connect()
        for attempt in range(0,self.tries):
            gz = ftp.ftp_open("/".join(path),  mode='rt', avoid_collision=True, delete=True)
            try:
                md5gz = rcf.md5(gz.name)
            except:
                raise IOError(f'''Could not open GBFF for {accession}.''')
            if md5 == md5gz:
                gz.assembly = accession
                break

        # Return file object
        return gz

    def genome_path(self, accession, assembly_reports=None):
        """
        Fetch the path of a genome at the NCBI FTP site.

        Usage:
          from rotifer.db.ncbi import ftp
          gc = ftp.GenomeCursor()
          path = gc.genome_path("GCA_900547725.1")
          print("/".join(path))

        Parameters:
          accession : string
            Genome accession
          assembly_reports: pandas DataFrame
            (Optional) NCBI's ASSEMBLY_REPORTS tables
            This dataframe can be loade using:

              from rotifer.db.ncbi as ncbi
              ar = ncbi.assemblies()
        
        Returns:
          A tuple of two strings, empty when the genome is not found
        """
        from rotifer.db.ncbi import NcbiConfig
        from rotifer.db.ncbi import ftp as ncbiftp
        ftp = ncbiftp.connection(tries=self.tries, timeout=self.timeout, cache=self.cache)
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
            path = ftp.ftp_ls(path)
            if path.empty:
                return ()
            path = path.query(f'name.str.contains("{accession}")')
            path = path.sort_values(['name'], ascending=False).iloc[0]
            path = path.target + "/" + path['name']

            # Retrieve GBFF path
            if not path:
                return ()
            path = ftp.ftp_ls(path)
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
          from rotifer.db.ncbi import ftp
          gc = ftp.GenomeCursor()
          contigs, assembly = gc.genome_report("GCA_900547725.1")

        Parameters:
          accession : string
            Genome accession
        
        Returns:
          A tuple containing a Pandas DataFrame and a Pandas Series
          
          The Pandas DataFrame lists the assembly's contigs

          The Series contains the assembly properties and is 
          similar to the table from rotifer.db.ncbi.assemblies()

        """
        from rotifer.db.ncbi import ftp as ncbiftp
        ftp = ncbiftp.connection(tries=self.tries, timeout=self.timeout, cache=self.cache)

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
        if len(path):
            report = ftp.ftp_ls(path[0])
        else:
            return ([],pd.DataFrame(columns=ar[0]))
        report = report[report.name.str.contains("_assembly_report.txt")]
        report = path[0] + "/" + report.name.iloc[0]
        report = ftp.ftp_open(report)

        inar = True
        for row in report:
            row = row.strip()
            if row == "#" or row[0:2] == "##":
                inar = False
            elif row[0:15] == "# Sequence-Name":
                sc.append(['assembly'] + row[2:].split("\t"))
            elif inar and row[0:2] == "# ":
                ar.append(row[2:].split(":", maxsplit=1))
            elif row[0] != '#':
                sc.append([accession] + row.split("\t"))

        ar = pd.DataFrame(ar[1:], columns=ar[0])
        ar.value = ar.value.str.strip()
        ar.column = ar.column.replace(arcolumn)
        sc = pd.DataFrame(sc[1:], columns=sc[0])
        ar = ar.set_index("column")

        return sc, ar

class GenomeFeaturesCursor(rotifer.db.methods.GenomeFeaturesCursor, GenomeCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes

    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> from rotifer.db.ncbi import ftp
    >>> gfc = ftp.GenomeFeaturesCursor(g)
    >>> df = gfc.fetchall()

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
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            timeout=10,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=tries, batch_size=batch_size, threads=threads, timeout=timeout, cache=cache, *args, **kwargs)
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, rotifer.db.parallel.GeneNeighborhoodCursor, GenomeFeaturesCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
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
      If set to True, neighborhood data for eukaryotic genomes
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
    timeout: integer
      Maximum amount of time, in seconds, to wait 
      for a connection to the server
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
            eukaryotes=False,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=True,
            tries=3,
            batch_size=config['batch_size'],
            threads = config["threads"] or _defaults['threads'],
            timeout=10,
            cache=rotifer.config['cache'],
            *args, **kwargs
        ):

        threads = threads or _defaults['threads']

        super().__init__(
            column = column,
            before = before,
            after = after,
            min_block_distance = min_block_distance,
            strand = strand,
            fttype = fttype,
            eukaryotes = eukaryotes,
            exclude_type = exclude_type,
            autopid = autopid,
            codontable = codontable,
            progress = progress,
            tries = tries,
            batch_size = batch_size,
            threads = threads,
            *args, **kwargs
        )
        self.timeout = timeout
        self.cache = cache

# Is this library being used as a script?
if __name__ == '__main__':
    pass
