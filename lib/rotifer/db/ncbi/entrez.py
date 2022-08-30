import sys
import rotifer
import pandas as pd
from Bio import Entrez, SeqIO
from rotifer.db.ncbi import NcbiConfig
Entrez.email = NcbiConfig["email"] if "email" in NcbiConfig else None
Entrez.api_key = NcbiConfig["api_key"] if "api_key" in NcbiConfig else None
logger = rotifer.logging.getLogger(__name__)

class cursor():
    """
    Fetch data from NCBI using accession numbers and EUtilities.
    """
    def __init__(self, database="protein", rettype="gbwithparts", retmode="text", tries=3, sleep_between_tries=1, batch_size=200):
        """
        Query and retrieve data from NCBI using EUtilities.

        Parameters
        ----------
        database: string, default 'protein'
            Valid NCBI sequence database
            Accpted values are protein and nucleotide
        rettype: string, default ```gbwithparts```
            Set the format for data returned by the NCBI.
        retmode: string, default 'text'
            Choose between raw text formats ("text") or XML text ("xml")
        tries: int, default 3
            Number of times to try downloading data
        sleep_between_tries: int, default 1
            Number of seconds to wait before tring a new download

        Returns
        -------

        Usage
        -----
        >>> from rotifer.db.ncbi import entrez
        >>> eutils = entrez.cursor(database="protein")
        >>> s = eutils.fetch("YP_009724395.1")
        """
        self.database = database
        self.rettype = rettype
        self.retmode = retmode
        self.tries = tries
        self.sleep_between_tries = sleep_between_tries
        self.batch_size = batch_size
        self.missing = set()

    def efetch(self, accessions):
        """
        Fetch sequences using EFetch.

        Usage
        -----
        from rotifer.db.ncbi import entrez
        eutils = entrez.cursor(database="protein", tries=3)
        p = eutils.fetch("YP_009724395.1")

        Parameters
        ----------
        accessions: (list of) strings
          NCBI accession numbers
        """
        # Reset query
        if not (isinstance(accessions,list) or isinstance(accessions,pd.Series)):
            accessions = [accessions]

        # Select format
        retformat = self.rettype
        if retformat in ["gbwithparts","gb","gp"]:
            retformat = "genbank"

        # Fetch
        batches = list(set(accessions))
        batches = [ batches[x:x+self.batch_size] for x in range(0,len(batches),self.batch_size) ]
        for batch in batches:
            yield Entrez.efetch(
                db = self.database,
                rettype = self.rettype,
                retmode = self.retmode,
                max_tries = self.tries,
                sleep_between_tries = self.sleep_between_tries,
                id = ",".join(batch),
            )

    def elink(self, accessions, to="nuccore", linkname=None):
        """
        Find related database entries via NCBI's EUtilities.

        Usage:
          # download from NCBI's FTP site
          from rotifer.db.ncbi import entrez
          eutils = entrez.cursor()
          a = eutils.elink("YP_009724395.1")

        Returns:
          Pandas DataFrame

        Parameters:
          accessions : string or list of strings
            NCBI accessions to search links for
          to: string
            Name of the target database
          linkname: string
            Type of link between dbfrom and dbto
            If not set, {dbfrom}_{dbto} is used
        """
        import pandas as pd
        from Bio import Entrez

        # Fix input
        if not isinstance(accessions,list):
            accessions = [accessions]
        if not linkname:
            linkname = self.database + "_" + to

        data = []
        for acc in accessions:
            try:
                raw = list(Entrez.read(Entrez.elink(dbfrom=self.database, linkname=linkname, id=acc)))
            except:
                logger.info(f'Entrez.elink failed for accession {acc}, dbfrom: {self.database}, dbto: {to}. Error: '+str(sys.exc_info()[0]))
                continue
            for d in raw:
                for x in d["LinkSetDb"]:
                    for y in x["Link"]:
                        data.append([acc, d["IdList"][0], d["DbFrom"], x["LinkName"], to, y["Id"]])
        data = pd.DataFrame(data, columns=["qacc", "quid", "dbfrom", "linkname", "dbto", "tuid"])

        return data
