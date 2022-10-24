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
from copy import deepcopy

import rotifer
import rotifer.db.core
from rotifer import GlobalConfig
from rotifer.db.ncbi import NcbiConfig
from rotifer.db.ncbi import utils as rdnu
import rotifer.db.ncbi.ftp as ncbiftp
from rotifer.genome.utils import seqrecords_to_dataframe
logger = rotifer.logging.getLogger(__name__)
DefaultBasePath = os.environ["ROTIFER_DATA"] if 'ROTIFER_DATA' in os.environ else "/databases"
DefaultBasePath = os.path.join(DefaultBasePath,"genomes")

# Classes

#class GenomeCursor(BaseCursor):
class GenomeCursor(ncbiftp.GenomeCursor):
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
    basepath: string
      Path to a mirror of the genomes section of the 
      NCBI FTP site. Contents are expected to be the
      same or a subset of the genomes directory.
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
            progress=False,
            tries=1,
            batch_size=None,
            threads=15,
            timeout=10,
            basepath = DefaultBasePath,
        ):
        super().__init__(
            progress=progress,
            tries=tries,
            batch_size=batch_size,
            threads=threads,
        )
        self.timeout = timeout
        self.basepath = basepath

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

        # find genome and download
        path = self.genome_path(accession, assembly_reports=assembly_reports)
        if len(path) == 0:
           return None

        # Download genome
        path = os.path.join(*path)
        gz = rcf.open_compressed(path,  mode='rt')

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
            path = [ path[i : i + 3] for i in range(0, len(path), 3) ]
            path = os.path.join(self.basepath,'all',*path)
            try:
                ls = os.listdir(path)
            except:
                raise FileNotFoundError(f'Directory {path} not found')
            ls = [ x for x in sorted(ls) if accession in x ]
            if len(ls):
                ls = ls[-1] # Expected to be the latest version of the target genome
            else:
                raise FileNotFoundError(f'No subdirectory for genome {accession} in diretocy {path}')
            path = os.path.join(path,ls)

            # Retrieve GBFF path
            ls = os.listdir(path)
            ls = [ x for x in sorted(ls) if '.gbff.gz' in x ]
            if len(ls):
                ls = ls[0] # Only one GBFF is expected
            else:
                raise FileNotFoundError(f'No GBFF for {àccession} in {path}')
            path = (path, ls[0])

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
        ar = [['column','value'], ['assembly', accession]]
        sc = []

        # Find report 
        path = self.genome_path(accession)
        if len(path):
            report = os.listdir(path[0])
        else:
            return ([],pd.DataFrame(columns=ar[0]))
        report = [ x for x in report if "_assembly_report.txt" in x ]
        if len(report):
            report = os.path.join(path[0],report[0])
        else:
            return ([],pd.DataFrame(columns=ar[0]))

        # Parse report
        report = open(report)
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
    basepath: string
      Path to a mirror of the genomes section of the 
      NCBI FTP site. Contents are expected to be the
      same or a subset of the genomes directory.
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
            basepath = DefaultBasePath,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
        ):
        super().__init__(
            basepath=basepath,
            progress=progress,
            tries=tries,
            batch_size=batch_size,
            threads=threads,
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

class GeneNeighborhoodCursor(GenomeFeaturesCursor, ncbiftp.GenomeFeaturesCursor):
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
    basepath: string
      Path to a mirror of the genomes section of the 
      NCBI FTP site. Contents are expected to be the
      same or a subset of the genomes directory.
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
            basepath = DefaultBasePath,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size=None,
            threads=15,
        ):
        super().__init__(
            basepath = basepath,
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
                    error = f'{sys.exc_info()[1]}'
                    if logger.getEffectiveLevel() <= logging.DEBUG:
                        logger.exception(error)
                    continue
                except ValueError:
                    error = f'{sys.exc_info()[1]}'
                    if logger.getEffectiveLevel() <= logging.DEBUG:
                        logger.exception(error)
                    break
                except:
                    error = f'{sys.exc_info()[1]}'
                    #error = f'Failed to download genome {accession}: {sys.exc_info()[1]}'
                    if logger.getEffectiveLevel() <= logging.DEBUG:
                        logger.exception(error)
                    continue

                if isinstance(stream, types.NoneType):
                    self._add_to_missing(expected, accession, error)
                    continue

                # Use parser to process results
                try:
                    obj = self.parser(stream, accession, assemblies[accession])
                    break
                except:
                    error = f"Failed to parse genome {accession}:"
                    if logger.getEffectiveLevel() <= logging.DEBUG:
                        logger.exception(error)

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
        if not self.eukaryotes:
            from rotifer.db.ncbi import entrez
            contigs, report = self.genome_report(accession)
            if len(report) > 0:
                tc = entrez.TaxonomyCursor()
                taxonomy = tc[report.loc['taxid'][0]]
                if taxonomy.loc[0,"superkingdom"] == "Eukaryota":
                    raise ValueError(f"Eukaryotic genome {accession} ignored.")
        stream = self.open_genome(accession)
        if isinstance(stream,types.NoneType):
            raise ValueError(f'Unable to access files for genome {accession}.')
        return stream

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
        # Make sure some content, even if empty, is always returned
        if len(result) == 0:
            result = [ seqrecords_to_dataframe([]) ]
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

    def fetchone(self, accessions, ipgs=None):
        """
        Asynchronously fetch gene neighborhoods from NCBI.

        Parameters
        ----------
        accessions: list of strings
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
        targets = deepcopy(accessions)
        if not isinstance(targets,typing.Iterable) or isinstance(targets,str):
            targets = [targets]
        targets = set(targets)
        todo = deepcopy(targets)

        # Make sure we have usable IPGs
        if isinstance(ipgs,types.NoneType):
            from rotifer.db.ncbi import entrez
            if self.progress:
                logger.warn(f'Downloading IPGs for {len(todo)} proteins...')
            size = self.batch_size
            ic = entrez.IPGCursor(progress=self.progress, tries=self.tries, threads=self.threads)
            ipgs = ic.fetchall(todo)
            #if len(ic.missing):
            #    self._add_to_missing(ic.missing.index.to_list(), np.nan, "No IPGs at NCBI")
        ipgids = set(ipgs[ipgs.pid.isin(todo) | ipgs.representative.isin(todo)].id)
        ipgs = ipgs[ipgs.id.isin(ipgids) & (ipgs.assembly.notna() | ipgs.nucleotide.notna())]

        # Check for proteins without IPGs
        missing = set(todo) - set(ipgs.pid).union(ipgs.representative)
        if missing:
            self._add_to_missing(missing,np.NaN,"Not found in IPGs")
            todo -= missing
        if len(ipgs) == 0:
            return [seqrecords_to_dataframe([])]

        # Select best IPGs
        assemblies = rdnu.best_ipgs(ipgs)

        # Mark proteins without assembly
        missing = assemblies[assemblies.assembly.isna()]
        if len(missing):
            for idx, row in missing.iterrows():
                if row['pid'] in todo:
                    acc = row['pid']
                elif row['representative'] in todo:
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
        genomes = set(assemblies.assembly.unique())
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            if self.progress:
                pids = set(assemblies.pid).union(assemblies.representative)
                pids = len(pids.intersection(targets))
                logger.warn(f'Downloading {len(genomes)} genomes for {pids} proteins...')
                p = tqdm(total=len(genomes), initial=0)
            tasks = []
            for chunk in self.splitter(assemblies):
                tasks.append(executor.submit(self.worker, chunk))
            completed = set()
            for x in as_completed(tasks):
                data = x.result()
                for s in data['missing'].iterrows():
                    if s[0] in targets:
                        self.missing.loc[s[0]] = s[1]
                for obj in data['result']:
                    found = self.getids(obj)
                    done = genomes.intersection(found)
                    if self.progress and len(done) > 0:
                        p.update(len(done))
                    genomes = genomes - done
                    completed.update(self._pids(obj))
                    self.missing.drop(completed, axis=0, inplace=True, errors='ignore')
                    yield obj

    def fetchall(self, accessions, ipgs=None):
        """
        Fetch genomes.

        Parameters
        ----------
        accessions: list of strings
          NCBI protein identifiers

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
        """
        stack = []
        for df in self.fetchone(accessions, ipgs=ipgs):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])

# Is this library being used as a script?
if __name__ == '__main__':
    pass
