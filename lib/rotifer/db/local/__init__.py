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
import rotifer.devel.beta.sequence as rdbs
logger = rotifer.logging.getLogger(__name__)

# Defaults
_config = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
    **loadConfig(__name__.replace('rotifer.',':'))
}

class esl_sfetch:
    def __init__(self, query=None, database_path=_config["local_database_path"], batch_size=200, threads=5):
        self.query = query
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

    def fetch(self, query=None):
        """
        Fetch many sequences from the database.

        Parameters
        ----------
        query: (list of) strings
          One or more sequence idenntifiers
        batch_size: int, default 200
          Number of sequences to process in each esl-sfetch call
        threads: int, default 20
          Number of batches (esl-sfetch calls) to run in parallel

        Returns
        -------
        A list of Bio.SeqRecords
        """
        import rotifer.devel.alpha.gian_func as gian_func
        if query:
            self.query = query
            self.missing = set()
        seqrecords = []
        for batch in gian_func.chunks(self.query, self.batch_size):
            seqrecords.extend(self._fetch(batch))
        return seqrecords
    
    def _fetch(self, query):
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
        with tempfile.NamedTemporaryFile(mode="wt") as tmp:
            tmp.write("\n".join([ str(x) for x in query ])+"\n")
            tmp.flush()
            s = ["esl-sfetch", "-f", self.path, tmp.name]
            s = subprocess.run(s, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if s.stderr:
                d = []
                for acc in query:
                    t = self[acc]
                    if t:
                        d.append(t)
                s = d
            else:
                s = [ self._clean_description(x) for x in SeqIO.parse(StringIO(s.stdout),"fasta") ]
            return s
