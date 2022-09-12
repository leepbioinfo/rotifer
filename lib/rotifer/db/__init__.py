import os
import rotifer
logger = rotifer.logging.getLogger(__name__)

# FUNCTIONS

def proteins(query, methods=['esl_sfetch','entrez'], local_database_path=os.path.join(rotifer.GlobalConfig["data"],"fadb","nr","nr"), entrez_database="protein", batch_size=200, threads=5, tries=3):
    """
    Fetch proteins sequences from local or remote databases.
    
    Note
    ----
    This method is optimized for speed and does not provide
    access to sequence annotation data, because all data is
    fetched as FASTA formatted data streams.

    Parameters
    ----------
    query: list of string
      Sequence identifiers
    methods: list of strings
      List of keywords listing methods of data retrieval.
      The following methods are available:
      - esl_sfetch
        rotifer.db.local.EaselFastaCursor
      - entrez
        rotifer.db.ncbi.entrez.FastaCursor
    local_database_path: strin
      Path to a (possibly indexed) sequence file.
      In most cases the local file will be a FASTA file.
    entrez_database: string
      Name the Entrez sequence database
    batch_size: int, default 200
      Number of sequences to retrieve per iteration
    threads: integer
      Number of simultaneous threads to run while fetching data
    tries: integer
      Maximum number of attempts to download from remote databases

    Returns
    -------
      One or a list of Bio.SeqRecord objects
    """
    import rotifer.db.local as localDB
    from rotifer.db.ncbi import entrez

    result = []
    targets = set(query)
    for method in methods:
        if method == 'esl_sfetch':
            cursor = localDB.EaselFastaCursor(database_path=local_database_path, batch_size=batch_size, threads=threads)
        elif method == 'entrez':
            cursor = entrez.FastaCursor(database=entrez_database, batch_size=batch_size, threads=threads, tries=tries)
        for s in cursor.fetch_each(targets):
            result.append(s)
        targets = cursor.missing
    if targets:
        logger.warn(f'A total of {len(targets)} sequences could not be found: targets')

    return result
