import os
import logging
from ftplib import FTP
from rotifer.core import GlobalConfig

# Load NCBI assembly reports
def ftp_get(target, avoid_collision=False, outdir=GlobalConfig['cache']):
    '''
    Download a file from the NCBI's FTP site.

    Usage:
      import rotifer.db.ncbi as ncbi
      localpath = ncbi.ftpget('genomes/README.txt')

    Parameters:
      target : string
        URL or relative path of the file
      avoid_collision : boolean, default False
        Avoid name collision by adding random suffixes
      outdir : string
        Output directory

    Returns:
      Path to the downloaded file
    '''
    from tempfile import NamedTemporaryFile
    from rotifer.db.ncbi import NcbiConfig

    # Set log format
    logger = logging.getLogger('rotifer.db.ncbi.ftp')
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

    # As we are about to start downloading our next target, make sure
    # we are still connected
    try:
        ftp.pwd()
    except:
        ftp = FTP(NcbiConfig['ftpserver'])
        ftp.login()

    # To avoid problems with very long names, I change to the target
    # directory and, later, back to /
    p = target.replace('ftp://'+NcbiConfig['ftpserver']+"/","").split("/")
    for i in list(range(len(p)-1)):
        ftp.cwd(p[i])
    ftp.retrbinary("RETR " + p[-1], outfh.write)
    ftp.cwd("/")
    outfh.close()

    # Return pandas object
    return outfile

# Load NCBI assembly reports
def ftp_open(target,  mode='rt', avoid_collision=True, delete=True, cache=GlobalConfig['cache']):
    '''
    Open a file stored at the NCBI FTP site.

    Note:
      Compressed files will be uncompressed on the fly.

    Usage:
      # Open a genome's GBF file
      import rotifer.db.ncbi as ncbi
      path = "genomes/all/GCA/900/547/725/GCA_900547725.1_UMGS1014/"
      path += "GCA_900547725.1_UMGS1014_genomic.gbff.gz"
      fh = ncbi.ftpopen(path, mode="rt")

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
      cache : string
        Path to a local directory to use as a temporary cache

    Returns:
      A data stream (file handle) is returned, similar to when a local
      file is accessed by `open()`.
    '''
    from tempfile import _TemporaryFileWrapper
    outfile = ftp_get(target, outdir=cache, avoid_collision=avoid_collision)
    outfh = _TemporaryFileWrapper(_hook_compressed_text(outfile, mode), outfile, delete)

    # Return pandas object
    return outfh

# Is this library being used as a script?
if __name__ == '__main__':
    pass
def _hook_compressed_text(filename, mode='r', encoding='utf8'):
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
