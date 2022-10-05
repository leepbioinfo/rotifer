import os
import sys
import types
import socket
import typing
import logging
import numpy as np
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

    # Download files
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
            This dataframe can be loaded using:

              from rotifer.db.ncbi as ncbi
              ar = ncbi.assemblies()

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
              ar = ncbi.assemblies()
        
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
        if len(path):
            report = self.ftp_ls(path[0])
        else:
            return None
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
                ar.append(row[2:].split(":", maxsplit=1))
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
        self.cache = cache

    def getids(self,obj):
        if isinstance(obj,types.NoneType):
            return set()
        if isinstance(obj,list):
            return set([ x.assembly for x in obj ])
        else:
            return {obj.assembly}

    def fetcher(self, accession):
        import rotifer.db.ncbi.ftp as ncbiftp
        gfc = ncbiftp.cursor(tries=1, cache=self.cache)
        return gfc.open_genome(accession)

    def parser(self, stream, accession):
        from Bio import SeqIO
        stack = []
        for s in SeqIO.parse(stream,"genbank"):
            s.assembly = accession
            stack.append(s)
        stream.close()
        return stack

    def worker(self,accessions):
        stack = []
        for accession in accessions:
            objlist = self[accession]
            if isinstance(objlist,types.Nonetype):
                continue
            if len(objlist) != 0:
                stack.extend(objlist)
        return stack

class GenomeFeaturesCursor(GenomeCursor):
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
            threads=threads,
            cache=cache,
        )
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

    def getids(self,obj):
        if isinstance(obj,types.NoneType):
            return set()
        if isinstance(obj,list):
            return set([ x.assembly for x in obj ])
        else:
            return set(obj.assembly)

    def parser(self, stream, accession):
        from Bio import SeqIO
        data = SeqIO.parse(stream,"genbank")
        data = seqrecords_to_dataframe(
            data,
            exclude_type = self.exclude_type,
            autopid = self.autopid,
            assembly = accession,
            codontable = self.codontable,
        )
        stream.close()
        return data

    def worker(self, accessions):
        stack = []
        for accession in accessions:
            df = self[accession]
            if len(df) != 0:
                stack.append(df)
        return stack

    def fetchall(self, accessions):
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
        for df in self.fetchone(accessions):
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
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
            cache=GlobalConfig['cache'],
        ):
        super().__init__(
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
        self.eukaryotes = eukaryotes
        self.missing = pd.DataFrame(columns=["noipgs","eukaryote","assembly","error",'class'])

    def _pids(self, obj):
        columns = ['pid']
        if 'replaced' in obj.columns:
            columns.append('replaced')
        ids = obj.melt(id_vars=['nucleotide'], value_vars=columns, value_name='id', var_name='type')
        ids.drop('type', axis=1, inplace=True)
        ids.set_index('nucleotide', inplace=True)
        ids.drop_duplicates(inplace=True)
        return ids.id.tolist()

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
            ipgs = ic.fetchall(protein)
        ipgs = ipgs[ipgs.id.isin(ipgs[ipgs.pid.isin(protein) | ipgs.representative.isin(protein)].id)]
        best = rdnu.best_ipgs(ipgs)
        best = best[best.assembly.notna()]
        ipgs = ipgs[ipgs.assembly.isin(best.assembly)]
        missing = set(protein) - set(ipgs.pid) - set(ipgs.representative)
        if missing:
            self._add_to_missing(missing,np.NaN,"No IPGs")
            return objlist

        # Identify DNA data
        assemblies, nucleotides = rdnu.ipgs_to_dicts(ipgs)

        # Download and parse
        objlist = []
        for accession in assemblies.keys():
            expected = set([ y for x in assemblies[accession].items() for y in x ])

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
                    error = f'Failed to download genome {accession}: {sys.exc_info()[1]}'
                    logger.debug(error)
                    continue

                if isinstance(stream, types.NoneType):
                    self._add_to_missing(expected, accession, error)
                    continue

                # Use parser to process results
                try:
                    obj = self.parser(stream, accession, assemblies[accession])
                    break
                except:
                    error = f"Failed to parse genome {accession}: {str(sys.exc_info()[1])}"
                    logger.debug(error)

            if isinstance(obj, types.NoneType):
                self._add_to_missing(expected, accession, error)
            elif len(obj) == 0:
                error = f'No anchors in genome {accession}'
                self._add_to_missing(expected, accession, error)
            else:
                objlist.append(obj)

        # No data?
        if len(objlist) == 0:
            return seqrecords_to_dataframe([])

        # Concatenate and evaluate
        objlist = pd.concat(objlist, ignore_index=True)

        # Return data
        if len(objlist) > 0:
            self.missing.drop(self._pids(objlist), axis=0, inplace=True, errors='ignore')
        return objlist

    def fetcher(self, accession):
        import rotifer.db.ncbi.ftp as ncbiftp
        cursor = ncbiftp.cursor(tries=1, cache=self.cache)
        ok = True
        if not self.eukaryotes:
            from rotifer.db.ncbi import entrez
            contigs, report = cursor.genome_report(accession)
            tc = entrez.TaxonomyCursor()
            taxonomy = tc[report.loc['taxid'][0]]
            if taxonomy.loc[0,"superkingdom"] == "Eukaryota":
                raise ValueError(f"Eukaryotic genome {accession} ignored.")
        return cursor.open_genome(accession)

    def parser(self, stream, accession, proteins):
        data = super().parser(stream, accession)
        data = data.neighbors(
            data[self.column].isin(proteins.keys()),
            before = self.before,
            after = self.after,
            min_block_distance = self.min_block_distance,
            strand = self.strand,
            fttype = self.fttype,
        )
        data['replaced'] = data.pid.replace(proteins)
        return data

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
            size = max(int(ipgs.assembly.nunique()/self.threads),1)
        batch = []
        for x, y in ipgs.groupby('assembly'):
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
        ipgs : Pandas dataframe
          This parameter may be used to avoid downloading IPGs
          from NCBI. Example:

          >>> from rotifer.db.ncbi import entrez
          >>> from rotifer.db.ncbi import ftp
          >>> ic = ncbi.IPGCursor(batch_size=1)
          >>> gnc = ftp.GeneNeighborhoodCursor(progress=True)
          >>> i = ic.fetchall(['WP_063732599.1'])
          >>> n = gnc.fetchall(['WP_063732599.1'], ipgs=i)

        Returns
        -------
        A generator for rotifer.genome.data.NeighborhoodDF objects
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Make sure no identifiers are used twice
        if not isinstance(proteins,typing.Iterable) or isinstance(proteins,str):
            proteins = [proteins]
        proteins = set(proteins)

        # Make sure we have usable IPGs
        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            if self.progress:
                logger.info(f'Downloading IPGs for {len(proteins)} proteins...')
            size = self.batch_size
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries, threads=self.threads)
            ipgs = ic.fetchall(list(proteins))
            if len(ic.missing):
                self._add_to_missing(ic.missing.index.to_list(), np.nan, "No IPGs at NCBI")
        ipgs = ipgs[ipgs.pid.isin(proteins) | ipgs.representative.isin(proteins)]
        ipgs = ipgs[ipgs.assembly.notna() | ipgs.nucleotide.notna()]

        # Check for proteins without IPGs
        missing = set(proteins) - set(ipgs.pid).union(ipgs.representative)
        if missing:
            self._add_to_missing(missing,np.NaN,"Not found in IPGs")
        if len(ipgs) == 0:
            return [seqrecords_to_dataframe([])]

        # Select best IPGs
        assemblies = rdnu.best_ipgs(ipgs)

        # Mark proteins without assembly
        missing = assemblies[assemblies.assembly.isna()]
        if len(missing):
            for idx, row in missing.iterrows():
                if row['pid'] in proteins:
                    acc = row['pid']
                elif row['representative'] in proteins:
                    acc = row['representative']
                else:
                    continue
                if pd.isna(row['nucleotide']):
                    error = f"No nucleotide or assembly for protein {acc}"
                else:
                    error = f"Fetch protein {acc} from nucleotide {row['nucleotide']}"
                self._add_to_missing(acc,row['assembly'],error)

        # filter good IPGs for the best assemblies
        assemblies = assemblies[assemblies.assembly.notna()]
        assemblies = ipgs[ipgs.assembly.isin(assemblies.assembly)]

        # Split jobs and execute
        todo = set(assemblies.assembly.unique())
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                targets = set(assemblies.pid).union(assemblies.representative)
                targets = len(targets.intersection(proteins))
                logger.info(f'Downloading {len(todo)} genomes for {targets} proteins...')
                p = tqdm(total=len(todo), initial=0)
            tasks = []
            for chunk in self.splitter(assemblies):
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
                    completed.update(self._pids(obj))
                    self.missing.drop(completed, axis=0, inplace=True, errors='ignore')
                    yield obj

    def fetchall(self, proteins, ipgs=None):
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
        for df in self.fetchone(proteins, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

# Is this library being used as a script?
if __name__ == '__main__':
    pass
