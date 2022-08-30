import os
import sys
import logging
import pandas as pd
from ftplib import FTP
import rotifer
from rotifer import GlobalConfig
from rotifer.db.ncbi import NcbiConfig
logger = rotifer.logging.getLogger(__name__)

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
        self.connection = FTP(url, timeout=timeout)
        self.connection.login()

    def connect(self):
        """
        Connect or reconnect to server.
        """
        try:
            self.connection.sendcmd("NOOP")
        except:
            try:
                self.connection = FTP(NcbiConfig["ftpserver"], timeout=self.timeout)
                self.connection.login()
            except:
                self._error = 1
                return
        self._error = 0

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
        attempt = 0
        self.connect()
        while attempt < self.tries:
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
            attempt += 1
        if not md5:
            self._error = 1
            return None

        # Download genome
        gz = None
        md5gz = None
        attempt = 0
        self.connect()
        while md5 != md5gz and attempt < self.tries:
            gz = self.ftp_open("/".join(path),  mode='rt', avoid_collision=True, delete=True)
            if self._error:
                logger.error(f'''Could not open GBFF for {accession}.''')
            else:
                md5gz = rcf.md5(gz.name)
                if md5 == md5gz:
                    break
            attempt += 1

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
            path = self.ftp_ls(path)
            if path.empty:
                return ()
            path = path[path['name'].str.contains(".gbff.gz")]
            if path.empty:
                return ()
            path = (path.target.iloc[0], path['name'].iloc[0])

        return path

# Is this library being used as a script?
if __name__ == '__main__':
    pass
