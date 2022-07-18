def elink(query, dbfrom="protein", dbto="nuccore", linkname=None, verbose=False):
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
      verbose : boolean
        Whether to print warnings and error messages
    """
    import pandas as pd
    from Bio import Entrez

    # Fix input
    if not isinstance(query,list):
        query = [query]
    if not linkname:
        linkname = dbfrom + "_" + dbto

    data = []
    for acc in query:
        try:
            raw = list(Entrez.read(Entrez.elink(dbfrom=dbfrom, linkname=linkname, id=acc)))
        except:
            if (verbose):
                print(f'{__name__}: Entrez.elink failed for accession {acc}, dbfrom: {dbfrom}, dbto: {dbto}. Error: '+str(sys.exc_info()[0]), file=sys.stderr)
            continue
        for d in raw:
            for x in d["LinkSetDb"]:
                for y in x["Link"]:
                    data.append([acc, d["IdList"][0], d["DbFrom"], x["LinkName"], dbto, y["Id"]])
    data = pd.DataFrame(data, columns=["qacc", "quid", "dbfrom", "linkname", "dbto", "tuid"])

    return data

