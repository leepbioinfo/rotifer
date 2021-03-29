#!/usr/bin/env python3

# Python's standard libraries
import os
import sys
import logging

# External libraries
import Bio.SeqIO
import pandas as pd

# Rotifer libraries and data
from rotifer.core    import GlobalConfig
from rotifer.db.ncbi import NcbiConfig

# Load NCBI genome assemblies
def genomes(ncbi, parser=Bio.SeqIO.parse, parser_args=['genbank'],
        parser_kwargs={}, outdir=GlobalConfig['cache'],
        fetch=['ftp'], assembly_reports=None, mode='r',
        concat=True, tempfile=True, delete=True, *args, **kwargs):
    '''
    Load NCBI assembly reports from a local directory or FTP.

    Usage:
      from rotifer.db.ncbi import ncbi
      a = ncbi()
      a.submit(['GCF_900504695.1', 'GCF_004636045.1', 'GCF_902726645.1'])
      b = a.parse(method='genomes', strategy=['ftp'])

    Parameters:
      ncbi             : rotifer.db.ncbi object
      parser           : parser function (default: Bio.SeqIO.parse)
      parser_args      : list of arguments to send to the parser (default: 'genbank')
      parser_kwargs    : list of named arguments to be consumed by the parser
      outdir           : where to store downloaded files
      concat           : (boolean) if True, return a single data stream
      tempfile         : avoid collision by adding a random string to file names
      delete           : if True, downloaded files are deleted when closed
      fetch            : list of rotifer.db.ncbi.fetch methods
      assembly_reports : rotifer.db.ncbi.read.assembly_reports dataframe
                         If not given, assembly reports are downloaded

    Returns:
      Whatever object the parser creates
      The default parser (Bio.SeqIO.parse) return Bio.SeqRecord objects
    '''
    #from rotifer.db.ncbi import ncbi as ncbiClass
    ncbiClass = ncbi.__class__

    # Check query list: should not be empty
    queries = ncbi.submit()
    if not queries:
        print('The list of queries is empty! Set it with submit() first!', file=sys.stderr)
        return None

    # Apply each fetcher and parse
    io = []
    for strategy in fetch:
        if strategy == 'ftp':
            if not isinstance(assembly_reports,pd.DataFrame):
                assembly_reports = ncbi.read('assembly_reports', baseurl = f'ftp://{NcbiConfig["ftpserver"]}/genomes/ASSEMBLY_REPORTS')
            files = assembly_reports[assembly_reports.assembly.isin(queries)]
            found = list(files.assembly) 
            [ ncbi.missing(x, replace=False) for x in queries if not ((x in found) or (x in ncbi.missing())) ]
            files = list(files.ftp_path.apply(lambda x: x.replace(f'ftp://{NcbiConfig["ftpserver"]}','') + '/' + os.path.basename(x) + '_genomic.gbff.gz'))
            if files:
                ncbiObj = ncbiClass(query=files)
                try:
                    fh = ncbiObj.fetch(method='ftp', outdir=outdir, mode=mode, concat=concat, tempfile=tempfile, delete=delete)
                    if concat:
                        io.append(fh)
                    else:
                        for f in fh:
                            io.append(f)
                except:
                    print(f'There was an error downloading assemblies from NCBI with strategy {strategy}: {sys.exc_info()[1]}', file=sys.stderr)
                ncbi.missing(ncbiObj.missing(), replace=False)
                ncbi.files.extend(ncbiObj.files)
        elif strategy == 'entrez':
            pass

    # Flatten results from different strategies and return
    # --> note: currently unsupported! Add append support to fileinput
    if parser and mode:
        io = [ parser(x, *parser_args, **parser_kwargs) for x in io ]
    if concat:
        if len(io) == 1:
            io = io[0]
    return io

# Is this library being used as a script?
if __name__ == '__main__':
    pass
