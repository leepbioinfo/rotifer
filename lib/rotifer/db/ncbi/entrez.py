import rotifer
import pandas as pd
from Bio import Entrez, SeqIO
from rotifer.db.ncbi import NcbiConfig
Entrez.email = NcbiConfig["email"] if "email" in NcbiConfig else None
Entrez.api_key = NcbiConfig["api_key"] if "api_key" in NcbiConfig else None
logger = rotifer.logging.getLogger(__name__)

class cursor():
    """
    This object represents a query to NCBI's EUtilities.
    """
    def __init__(self, query=[], database="protein", rettype="gbwithparts", retmode="text", tries=3, sleep_between_tries=1, batch_size=200):
        """
        Query and retrieve data from NCBI using EUtilities.

        Usage
        -----
        from rotifer.db.ncbi import entrez
        eutils = entrez.cursor("YP_009724395.1", tries=3)
        for s in eutils:
            print(s.id)

        Parameters
        ----------
        query: (list of) strings
            NCBI sequence accession numbers
        database: string, default 'protein'
            Valid NCBI sequence database
            Accpted values are protein and nucleotide
        rettype: string, default ```text```
            Set the format of sequences returned by the NCBI.
            The NCBI database will return raw text to be
            parsed by Bio.SeqIO. Valid formats are fasta
            and gbwithparts (classic GenBank format).
        retmode: string, default 'text'
            Choose between classic raw text formats ("text"), such
            as fasta or genbank, or XML text ("xml").
        tries: int, default 3
            Number of times to try downloading data
        sleep_between_tries: int, default 1
            Number of seconds to wait before tring a new download
        """
        if not (isinstance(query,list) or isinstance(query,pd.Series)):
            query = [query]
        self.query = query
        self.database = database
        self.rettype = rettype
        self.retmode = retmode
        self.tries = tries
        self.sleep_between_tries = sleep_between_tries
        self.batch_size = batch_size
        self.missing = set()

    def fetch(self, query=None, rettype=None):
        """
        Fetch sequences using EFetch.

        Usage
        -----
        from rotifer.db.ncbi import entrez
        eutils = entrez.cursor("YP_009724395.1", tries=3)
        for s in eutils:
            print(s.id)

        Parameters
        ----------
        query: (list of) strings
            NCBI sequence accession numbers
            Note: use of this argument replaces the
                  internal list of sequence accessions.
        """
        # Reset query
        if query:
            if not (isinstance(query,list) or isinstance(query,pd.Series)):
                query = [query]
            self.query = query
            self.missing = set()

        # Select format
        retformat = self.rettype
        if retformat in ["gbwithparts","gb","gp"]:
            retformat = "genbank"

        # Fetch
        stack = []
        batches = [ self.query[x:x+self.batch_size] for x in range(0,len(self.query), self.batch_size) ]
        for batch in batches:
            io = Entrez.efetch(
                db = self.database,
                rettype = self.rettype,
                retmode = self.retmode,
                max_tries = self.tries,
                sleep_between_tries = self.sleep_between_tries,
                id = ",".join(batch),
            )
            io = SeqIO.parse(io,retformat)
            try:
                stack.extend([ s for s in io ])
            except:
                self.missing = self.missing.union(batch)
        return stack

    def elink(query=None, dbfrom="protein", dbto="nuccore", linkname=None):
        """
        Find related database entries via NCBI's EUtilities.

        Usage:
          # download from NCBI's FTP site
          import rotifer.db.ncbi as ncbi
          a = ncbi.elink("YP_009724395.1")

        Returns:
          Pandas DataFrame

        Parameters:
          query : string or list of strings
            NCBI accessions to search links for
          dbfrom : string
            Name of the input database
          dbto: string
            Name of the target database
          linkname: string
            Type of link between dbfrom and dbto
            If not set, {dbfrom}_{dbto} is used
        """
        import pandas as pd
        from Bio import Entrez

        # Fix input
        if not query:
            query = self.query
        if not isinstance(query,list):
            query = [query]
        if not linkname:
            linkname = dbfrom + "_" + dbto

        data = []
        for acc in query:
            try:
                raw = list(Entrez.read(Entrez.elink(dbfrom=dbfrom, linkname=linkname, id=acc)))
            except:
                logger.info(f'Entrez.elink failed for accession {acc}, dbfrom: {dbfrom}, dbto: {dbto}. Error: '+str(sys.exc_info()[0]))
                continue
            for d in raw:
                for x in d["LinkSetDb"]:
                    for y in x["Link"]:
                        data.append([acc, d["IdList"][0], d["DbFrom"], x["LinkName"], dbto, y["Id"]])
        data = pd.DataFrame(data, columns=["qacc", "quid", "dbfrom", "linkname", "dbto", "tuid"])

        return data
