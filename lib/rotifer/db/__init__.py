import os
import rotifer
from . import local as localDB
from .ncbi import entrez
logger = rotifer.logging.getLogger(__name__)

# FUNCTIONS

def fetch_proteins(query, local_database_path=os.path.join(rotifer.GlobalConfig["data"],"fadb","nr","nr"), remote_database="entrez", batch_size=200, threads=5, tries=3):
    """
    Fetch proteins from local or remote databases

    Parameters
    ----------
    query: list of string
      Sequence identifiers
    local_database_path: string
      Path to a (possibly indexed) sequence file.
      In most cases the local file will be a FASTA file.
    remote_database: string
      Name of a remote sequence database resource
    batch_size: int, default 200
      Number of sequences to retrieve per iteration
    threads: integer
      Number of simultaneous threads to run while fetching data
    tries: integer
      Maximum number of attempts to download from remote databases
    """
    missing = set()
    easel = localDB.esl_sfetch(query, database_path=local_database_path, batch_size=batch_size, threads=threads)
    found = easel.fetch()
    if easel.missing:
        remoteDB = entrez.cursor(list(easel.missing), tries=tries)
        found.extend(remoteDB.fetch())
        if remoteDB.missing:
            logger.warn(f'A total of {len(remoteDB.missing)} sequences could not be found.')

    return found
