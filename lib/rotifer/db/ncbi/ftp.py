import os
import sys
import logging
import pandas as pd
from ftplib import FTP
from rotifer.core import GlobalConfig
from rotifer.db.ncbi import NcbiConfig

class NcbiFTPCursor():
    """
    This class represents a live connection to the NCBI FTP server.
    """
    def __init__(self, url=NcbiConfig['ftpserver'], tries=3, timeout=60, loglevel=None):
        self.url = url
        self.tries = tries
        self.connection = FTP(url, timeout=timeout)
        self.connection.login()

        # Create logger
        logger = logging.getLogger(f'{__name__}')
        shandler = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s : %(message)s")
        shandler.setFormatter(formatter)
        logger.addHandler(shandler)
        if loglevel:
            if isinstance(loglevel,str):
                logger.setLevel(getattr(logging,loglevel))
            else:
                logger.setlevel(loglevel)
        else:
            logger.setLevel(logging.ERROR)
        self.logger = logger

    def connect(self):
        """
        Connect or reconnect to server.
        """
        try:
            self.connection.sendcmd("NOOP")
        except:
            self.connection = FTP(NcbiConfig["ftpserver"], timeout=self.timeout)
            self.connection.login()

    # Load NCBI assembly reports
    def ftp_get(self, target, avoid_collision=False, outdir=GlobalConfig['cache']):
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
        self.logger.debug(f'started download of {NcbiConfig["ftpserver"]}{target}')

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

        # To avoid problems with very long names, I change to the target
        # directory and, later, back to /
        p = target.replace('ftp://',"")
        p = p.replace(NcbiConfig['ftpserver'],"")
        p = p[1:] if p[0] == "/" else p
        p = p.split("/")
        self.connect()
        for i in list(range(len(p)-1)):
            self.connection.cwd(p[i])
        self.connection.retrbinary("RETR " + p[-1], outfh.write)
        self.connection.cwd("/")
        outfh.close()

        # Return pandas object
        self.logger.debug(f'downloaded {NcbiConfig["ftpserver"]}{target}')
        return outfile

    # List files in ftp directory
    def ftp_ls(self, targets):
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
                self.logger.debug(f'''Could not retrive list for directory {target} at the NCBI's FTP site.''')
                continue
        d = pd.DataFrame(d)
        return d

    # Load NCBI assembly reports
    def ftp_open(self, target,  mode='rt', avoid_collision=True, delete=True, cache=GlobalConfig['cache']):
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
        outfile = self.ftp_get(target, outdir=cache, avoid_collision=avoid_collision)
        outfh = _TemporaryFileWrapper(self._hook_compressed_text(outfile, mode), outfile, delete)

        # Return pandas object
        return outfh

    def open_genome(self, accession, assembly_reports=None, cache=GlobalConfig['cache']):
        """
        Open the GBFF file of a genome from the NCBI FTP site.
        """
        import rotifer.core.functions as rcf
        path = self._genome_path(accession, assembly_reports=assembly_reports)
        if len(path) > 0:
            md5 = "/".join([path[0],"md5checksums.txt"])
            md5 = self.ftp_open(md5, cache=cache)
            md5 = pd.read_csv(md5, sep=' +', names=['md5','filename'])
            md5 = md5[md5.filename.str.contains('_genomic.gbff.gz')].md5.iloc[0]
            gz = None
            gzmd5 = None
            attempt = 1
            while md5 != gzmd5 and attempt < self.tries:
                gz = self.ftp_open("/".join(path),  mode='rt', avoid_collision=True, delete=True, cache=cache)
                gzmd5 = rcf.md5(gz.name)
                if md5 == gzmd5:
                    break
                attempt += 1
            return gz
        else:
            return None

    def _genome_path(self, accession, assembly_reports=None):
        from rotifer.db.ncbi import NcbiConfig
        path = ()

        # Extract genome path from assembly reports
        if isinstance(assembly_reports, pd.DataFrame) and not assembly_reports.empty:
            path = assembly_reports.query(f'assembly == "{accession}"').ftp_path.iloc[0]
            path = path.replace(f'ftp://{NcbiConfig["ftpserver"]}','')
            path = (path,os.path.basename(path) + "_genomic.gbff.gz")

        # Retrieve genome path for newest version
        else:
            path = accession.replace("GCF_","").replace("GCA_","")
            path = path.split(".")[0]
            path = "/".join([ path[i : i + 3] for i in range(0, len(path), 3) ])
            path = f'/genomes/all/{accession[0:3]}/{path}'
            path = self.ftp_ls(path)
            if path.empty:
                return ()
            path = path.query('type == "dir"')
            path = path.sort_values(['name'], ascending=False).iloc[0]
            path = path.target + "/" + path['name']

            # Retrieve GBFF path
            if not path:
                return ()
            self.connect()
            path = self.ftp_ls(path)
            if path.empty:
                return ()
            path = path[path['name'].str.contains(".gbff.gz")]
            if path.empty:
                return ()
            path = (path.target.iloc[0], path['name'].iloc[0])

        return path

    def _hook_compressed_text(self, filename, mode='r', encoding='utf8'):
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
