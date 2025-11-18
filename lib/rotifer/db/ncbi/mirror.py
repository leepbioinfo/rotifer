import os
import sys
import types
import typing
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from copy import deepcopy

import rotifer
import rotifer.db.parallel
import rotifer.db.methods
from rotifer import GlobalConfig
from rotifer.db.ncbi import config as NcbiConfig
from rotifer.db.ncbi import utils as rdnu
from rotifer.core.functions import loadConfig
from rotifer.genome.utils import seqrecords_to_dataframe
logger = rotifer.logging.getLogger(__name__)

# Configuration
from rotifer.core.functions import loadConfig
_defaults = {
    "path": NcbiConfig['mirror'] or os.path.join(rotifer.config['data'],"genomes"),
    "batch_size": None,
    "threads": int(np.floor(os.cpu_count()/2)),
}
config = loadConfig(__name__, defaults = _defaults)

# Classes

class GenomeCursor(rotifer.db.methods.GenomeCursor, rotifer.db.parallel.SimpleParallelProcessCursor):
    """
    Fetch genome sequences from a local mirror of the 
    NCBI genomes repository.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
    >>> from rotifer.db.ncbi import mirror
    >>> gc = mirror.GenomeCursor()
    >>> genomes = gc.fetchall(accessions)

    Parameters
    ----------
    path: string
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
            batch_size = config["batch_size"],
            threads = config["threads"] or _defaults['threads'],
            path = config["path"],
            *args, **kwargs):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=1, batch_size=batch_size, threads=threads)
        self.path = path

    def open_genome(self, accession, assembly_reports=None):
        """
        Open the GBFF file of a genome from a NCBI local mirror.

        Usage:
          # Just open
          
          from rotifer.db.ncbi import mirror
          gc = mirror.GenomeCursor()
          fh = gc.open_genome("GCA_900547725.1")

          # Open and parse to a list of Bio.SeqRecord

          from rotifer.db.ncbi import mirror
          from Bio import SeqIO
          gc = mirror.GenomeCursor()
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
        gz.assembly = accession

        # Return file object
        return gz

    def genome_path(self, accession, assembly_reports=None):
        """
        Fetch the path of a genome in the local NCBI mirror.

        Usage:
          from rotifer.db.ncbi import mirror
          gc = mirror.GenomeCursor()
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
                path = path.replace(f'ftp://{NcbiConfig["ftpserver"]}/genomes/','')
                path = (os.path.join(self.path,path),os.path.basename(path) + "_genomic.gbff.gz")

        # Retrieve genome path for newest version
        if len(path) == 0:
            path = accession[0:accession.find(".")].replace("_","")
            path = [ path[i : i + 3] for i in range(0, len(path), 3) ]
            path = os.path.join(self.path,'all',*path)
            if os.path.exists(path):
                ls = os.listdir(path)
            else:
                raise FileNotFoundError(f'No directory {path} for {accession}')
            ls = [ x for x in sorted(ls) if accession in x and os.path.isdir(os.path.join(path,x)) ]
            if len(ls):
                ls = ls[-1] # Expected to be the latest version of the target genome
            else:
                raise FileNotFoundError(f'Empty directory for {accession} in {path}')
            path = os.path.join(path,ls)

            # Retrieve GBFF path
            try:
                ls = os.listdir(path)
            except:
                raise IOError(f'Unable to read directory for {accession} in {path}')
            ls = [ x for x in sorted(ls) if '.gbff.gz' in x ]
            if len(ls):
                ls = ls[0] # Only one GBFF is expected
            else:
                ls = None
            path = (path, ls)

        return path

    def genome_report(self, accession):
        """
        Fetch genome assembly reports.

        Usage:
          from rotifer.db.ncbi import mirror
          gc = mirror.GenomeCursor()
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

class GenomeFeaturesCursor(rotifer.db.methods.GenomeFeaturesCursor, GenomeCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes

    >>> g = ['GCA_018744545.1', 'GCA_901308185.1']
    >>> from rotifer.db.ncbi import mirror
    >>> gfc = mirror.GenomeFeaturesCursor(g)
    >>> df = gfc.fetchall()

    Parameters
    ----------
    path: string
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
            path = config["path"],
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=1,
            batch_size = config["batch_size"],
            threads = config["threads"] or _defaults['threads'],
            *args, **kwargs
        ):
        threads = threads or _defaults['threads']
        super().__init__(progress=progress, tries=1, batch_size=batch_size, threads=threads, path=path, *args, **kwargs)
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, rotifer.db.parallel.GeneNeighborhoodCursor, GenomeFeaturesCursor):
    """
    Fetch genome annotation as dataframes.

    Usage
    -----
    Load a random sample of genomes, except eukaryotes
    >>> from rotifer.db.ncbi import mirror
    >>> gfc = mirror.GeneNeighborhoodCursor(progress=True)
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
    path: string
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
            path = config["path"],
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            tries=3,
            batch_size = config["batch_size"],
            threads = config["threads"] or _defaults['threads'],
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
        self.path = path

# Is this library being used as a script?
if __name__ == '__main__':
    pass
