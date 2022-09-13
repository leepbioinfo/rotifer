__doc__ = """
Rotifer connections to local databases
======================================
"""

import re
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from io import StringIO

import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
import rotifer.devel.beta.sequence as rdbs
logger = rotifer.logging.getLogger(__name__)

# Defaults
_config = loadConfig(__name__.replace('rotifer.',':'))
if not _config:
    _config = {}
_config = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
    **_config
}

class EaselFastaCursor:
    """
    Fetch biomolecular sequences using Easel's esl-sfetch.

    Parameters
    ----------
    database_path: string
      Path to a FASTA file.
    batch_size: int
      Number of sequence identifiers to process per thread
    threads: int
      Number of threads to run simultaneously

    """
    def __init__(self, database_path=_config["local_database_path"], batch_size=200, threads=5):
        self.path = database_path
        self.batch_size = batch_size
        self.threads = threads
        self.missing = set()
        if not os.path.exists(self.path):
            logger.error(f'{database_path}: no such file!')
            return None
        if not os.path.exists(self.path + ".ssi"):
            logger.warn(f'Building esl_sfetch index for {database_path}...')
            subprocess.run(["esl-sfetch","--index",self.path])

    def _clean_description(self, seqrec):
        seqrec.description = re.sub("\x01.+", "", seqrec.description.replace(seqrec.id, "").lstrip())
        return seqrec

    def __getitem__(self, accession):
        """
        Fetch one sequence from the database.

        Returns
        -------
        A Bio.SeqRecord object
        """
        sequence = subprocess.run(["esl-sfetch",self.path,accession], capture_output=True)
        if sequence.stderr:
            logger.warn(f'Esl-sfetch failed: no seq {accession} in database {self.path}')
            self.missing.add(accession)
        else:
            sequence = sequence.stdout.decode()
            sequence = SeqIO.read(StringIO(sequence),"fasta")
            sequence = self._clean_description(sequence)
            return sequence

    def fetch_each(self, query):
        for x in self.fetch_all(query):
            yield x

    def fetch_all(self, query):
        """
        Fetch many sequences from the database.

        Parameters
        ----------
        query: (list of) strings
          One or more sequence idenntifiers

        Returns
        -------
        A list of Bio.SeqRecords
        """
        query = set(query)
        import rotifer.devel.alpha.gian_func as gian_func
        seqrecords = []
        for batch in gian_func.chunks(list(query), self.batch_size):
            seqrecords.extend(self._fetch_many(batch))
        return seqrecords

    def _fetch_many(self, query):
        """
        Fetch many sequences from the database.

        Parameters
        ----------
        query: (list of) strings
          One or more sequence idenntifiers

        Returns
        -------
        A list of Bio.SeqRecords
        """
        import tempfile
        import subprocess
        from Bio import SeqIO
        from io import StringIO

        result = []
        query = set(query)
        with tempfile.NamedTemporaryFile(mode="wt") as tmp:
            tmp.write("\n".join([ str(x) for x in query ])+"\n")
            tmp.flush()
            s = ["esl-sfetch", "-f", self.path, tmp.name]
            s = subprocess.run(s, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if s.stderr:
                for acc in query:
                    t = self[acc]
                    if t:
                        result.append(t)
            else:
                found = set()
                for x in SeqIO.parse(StringIO(s.stdout),"fasta"):
                    x = self._clean_description(x)
                    found.add(x.id)
                    result.append(x)
                self.missing = self.missing.union(query - found)
        return result
