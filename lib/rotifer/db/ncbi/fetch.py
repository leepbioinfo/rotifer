#!/usr/bin/env python3

import os
import logging
from tempfile import NamedTemporaryFile, _TemporaryFileWrapper
from ftplib import FTP
from rotifer.io import fileinput
from rotifer.core import GlobalConfig
from rotifer.db.ncbi import NcbiConfig

# Load NCBI assembly reports
def ftp(ncbi, outdir=GlobalConfig['cache'], mode='r', concat=False, tempfile=False, delete=False, *args, **kwargs):
    '''
    Download data from NCBI's FTP site

    Usage:
      a = ncbi(query=['genomes/README.txt'])
      path = a.fetch(method="ftp")
      fh = a.fetch(method="ftp", mode="rt")

    Parameters:
      ncbi     : rotifer.db.ncbi.ncbi object
      outdir   : output directory
      mode     : set read mode ('r', 'rt', 'rb', etc) to open files

                 Note: if mode is set to False, files will not be 
                 open and, instead of file objects, fully qualified
                 paths to downloaded files will be returned.

      concat   : return one rotifer.io.fileinput.FileInput object
                 
                 This option automatically decompresses gzip and
                 bzip2 files.
                 
                 If True, implies mode='r'.

                 If set to False, a list of files (mode=False) or
                 file objects (mode='r') is returned

      tempfile : avoid name collision with temporary files
      delete   : files will be automatically removed when closed

    Returns:
      One or more file names or file objects.
      File objects may be file or rotifer.io.fileinput objects.
    '''

    # Set log format
    logger = logging.getLogger('rotifer.db.ncbi.fetch')
    logger.setLevel(logging.DEBUG)
    logger.info(f'main: downloading files from {NcbiConfig["ftpserver"]}')

    # Connect if necessary
    ftp = FTP(NcbiConfig['ftpserver'])
    ftp.login()

    # Create output directory, if necessary
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except:
            print(f'Error while trying to create download directory {outdir}')
            return None

    # Retrieve contents for each folder
    files = []
    for target in ncbi.submit():
        # Download
        if tempfile:
            parts = os.path.splitext(os.path.basename(target))
            prefix = parts[0] + '.' 
            suffix = None if parts[1] == '' else parts[1]
            outfh = NamedTemporaryFile(mode='w+b', suffix=suffix, prefix=prefix, dir=outdir, delete=False)
            outfile = outfh.name
        else:
            outfile = os.path.join(outdir, os.path.basename(target))
            outfh   = open(outfile,'wb')
        ftp.retrbinary("RETR " + target, outfh.write)
        outfh.close()

        # Autoopen (mode=True) and/or concatenate (concat=True)
        if concat:
            files.append(outfile)
        else:
            outfh = _TemporaryFileWrapper(hook_compressed_text(outfile, mode), outfile, delete)
            if mode:
                ncbi.files.extend([outfile])
                files.append(outfh)
            else:
                ncbi.files.extend([outfh])
                files.append(outfile)

    # Concatenate?
    if concat:
        ncbi.files.extend(files)
        files = fileinput.input(files, mode=mode, openhook=hook_compressed_text, delete=delete)

    # Return pandas object
    return files

def genomes(ncbi, outdir=GlobalConfig['cache'], assembly_reports=None, mode=False, concat=False, tempfile=False, delete=False, *args, **kwargs):
    '''
    Download NCBI genomes from FTP site to local directory.

    Usage:
      a = ncbi(query=['GCF_900504695.1', 'GCF_004636045.1', 'GCF_902726645.1'])
      p = a.fetch(method="genomes")

    Returns:
      One, if concat=True, or a list of file objects or file paths,
      if mode=False

    Parameters:
      ncbi             : rotifer.db.ncbi.ncbi object
      outdir           : output directory
      assembly_reports : rotifer.db.ncbi.read.assembly_reports dataframe
      mode             : set read mode ('r','rt','rb',etc) to open files

        Note: if mode is set to False, files will not be 
        open and, instead of file objects, fully qualified
        paths to downloaded files will be returned.

      concat   : return one rotifer.io.fileinput.FileInput object
                 
                 This option automatically decompresses gzip and
                 bzip2 files.
                 
                 If True, implies mode='r'.

                 If set to False, a list of files (mode=False) or
                 file objects (mode='r') is returned

      tempfile : avoid name collision with temporary files
      delete   : files will be automatically removed when closed
    '''
    return ncbi.parse('genomes', assembly_reports=assembly_reports, outdir=outdir, parser=None, mode=mode, concat=concat, tempfile=tempfile, delete=delete, *args, **kwargs)

# Internal method to use when returning concatenated file
# streams (handlers) with fileinput
def hook_compressed_text(filename, mode='r', encoding='utf8'):
    ext = os.path.splitext(filename)[1]
    if (mode == 'r') or not mode:
        mode = 'rt'
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode, encoding=encoding)
    elif ext == '.bz2':
        import bz2
        return bz2.open(filename, mode, encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)

if __name__ == '__main__':
    pass
