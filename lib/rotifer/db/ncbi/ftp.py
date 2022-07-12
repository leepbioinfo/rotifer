import os
import logging
import pandas as pd
from ftplib import FTP
from rotifer.core import GlobalConfig

# Load NCBI assembly reports
def ftp_get(target, avoid_collision=False, outdir=GlobalConfig['cache']):
    '''
    Download a file from the NCBI's FTP site.

    Usage:
      import rotifer.db.ncbi as ncbi
      localpath = ncbi.ftp_get('genomes/README.txt')

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

# List files in ftp directory
def ftp_ls(targets):
    '''
    List contents of an FTP directory.

    Usage:
      from rotifer.db.ncbi import ftp
      contents = ftp.ftp_ls('genomes')

    Parameters:
      target : (list of) string(s)
        Path of one or more diretory(ies)

    Returns:
      Pandas DataFrame
    '''
    from rotifer.db.ncbi import NcbiConfig
    import pandas as pd

    # Connect if necessary
    ftp = FTP(NcbiConfig['ftpserver'])
    ftp.login()

    # Process targets
    d = []
    if not (isinstance(targets,list) or isinstance(targets,tuple)):
        targets = [targets]
    for target in targets:
        try:
            for x in ftp.mlsd(target):
                if x[0] == "." or x[0] == "..":
                    continue
                x[1]["target"] = target
                x[1]["name"] = x[0]
                d.append(x[1])
        except:
            print(f'''Could not retrive list for directory {target} at the NCBI's FTP site.''')
            continue
    d = pd.DataFrame(d)
    return d

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

def open_genome(accession, assembly_reports=None, cache=GlobalConfig['cache']):
    """
    Open the GBFF file of a genome hosted at NCBI's FTP site.
    """
    path = _genome_path(accession, assembly_reports=assembly_reports)
    if path:
        return ftp_open(path,  mode='rt', avoid_collision=True, delete=True, cache=cache)
    else:
        return None

def _genome_path(accession, assembly_reports=None):
    from rotifer.db.ncbi import NcbiConfig
    path = None

    # Extract genome path from assembly reports
    if isinstance(assembly_reports, pd.DataFrame) and not assembly_reports.empty:
        path = assembly_reports.query(f'assembly == "{accession}"').ftp_path.iloc[0]
        path = path.replace(f'ftp://{NcbiConfig["ftpserver"]}','')

    # Retrieve genome path for newest version
    else:
        path = accession.replace("GCF_","").replace("GCA_","")
        path = path.split(".")[0]
        path = "/".join([ path[i : i + 3] for i in range(0, len(path), 3) ])
        path = f'/genomes/all/{accession[0:3]}/{path}'
        path = ftp_ls(path)
        if path.empty:
            return None
        path = path.query('type == "dir"')
        path = path.sort_values(['modify'], ascending=False).iloc[0]
        path = path.target + "/" + path['name']

    # Retrieve GBFF path
    if not path:
        return None
    path = ftp_ls(path)
    if path.empty:
        return None
    path = path[path['name'].str.contains(".gbff.gz")]
    if path.empty:
        return None
    path = (path.target + "/" + path['name']).iloc[0]
    return path

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

# Is this library being used as a script?
if __name__ == '__main__':
    pass
