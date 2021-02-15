#!/usr/bin/env python3

# Core libraries
import itertools
import os
import re
import sys
import types

# Pandas
import numpy as np
import pandas as pd

# Biopython
from Bio import SeqIO

def _get_basename(path):
    # Parse filename
    bn = path
    if path == sys.stdin:
        bn = "__stdin__"
    elif os.path.exists(path):
        bn = os.path.splitext(os.path.basename(path))[0]
    else:
        raise TypeError("Unknown file-like object")
    return bn

def SequenceNameToAssembly(path, s):
    return s.name

def GenbankAssemblyFromFilename(path, s):
    """
    Parse filenames with the regular expression
    
    ^(GC[AF]_\d+\.\d+)_\S+\.gbff(\.gz)?
    
    to extract Genbank assembly IDs.
    """
    aid = _get_basename(path)
    if re.match(r'^(GC[AF]_\d+\.\d+)',aid):
        aid = re.sub('^(GC[AF]_\d+\.\d+).*','\\1',aid)
    return aid

def parse(handle, informat, *args, assembly=GenbankAssemblyFromFilename, **kwargs):
    """
    Load genomic data from sequence file(s).
    This method crates an iterator capable of return SeqRecords.

    Arguments:
     - handle   : one or more file handles or filenames
     - informat : either gff or a Bio.SeqIO supported format

     - assembly : rule for setting genome assembly accession ID
                  Assembly identifiers are not part of the INSDC 
                  (http://www.insdc.org/) standard but can be
                  found in some NCBI entries under the DBLINK
                  tag, such as in many RefSeq entries:

                  DBLINK    BioProject: PRJNA224116
                            BioSample: SAMEA3138382
                            Assembly: GCF_000210855.2 <-- here!

                  This parameter allows the user to automate the definition
                  of genome identifiers. You can use either a string or a
                  lambda function that takes two arguments: file path and
                  sequence record object.

                  See GenbankAssemblyFromFilename and SequenceNameToAssembly
                  for examples of function compatible with this parameter.

     - retvalue : object(s) type to return. Supported objects are

     Other arguments depend on the parser (see methods below)
    """
    if not isinstance(handle,list):
        handle = [handle]

    # Process each file, return and wait for next request
    for eachhandle in handle:
        # Choose parser
        if informat == "gff":
            generator = gff(eachhandle, *args, **kwargs)
        else:
            generator = SeqIO.parse(eachhandle,informat)

        # Parse file basename
        bn = _get_basename(eachhandle)

        # Process each Bio.SeqRecord
        for s in generator:
            # Change number-only sequence names
            if re.match('^\d+$',s.id):
                s.id = bn + ":" + s.id
                s.name = s.id

            # Fix DESCRIPTION
            if not s.description:
                s.description = s.name + ", whole genome shotgun sequence"

            # Fix assembly
            if not any([ "Assembly:" == x[0:9] for x in s.dbxrefs ]):
                if assembly:
                    if isinstance(assembly,types.FunctionType):
                        assemblyStr = str(assembly(eachhandle,s))
                    else:
                        assemblyStr = assembly
                    s.dbxrefs.append(f'Assembly:{assemblyStr}')
                else:
                    s.dbxrefs.append(f'Assembly:{s.name}')

            # Return and wait next call
            yield s

def gff(path, seqid=None, organism=None, strain=None, description=None, topology='linear', keywords=[], ete3=None):
    """
    Parse GFF file

    Arguments:
     - path        : file path or file handle
     - seqid       : rename sequence identifiers using a user-defined function
                     This function will be given two arguments:
                      - the sequence object (a Bio.SeqRecord.SeqRecord object)
                      - the input file path
     - description : set sequence description
                     Could be a lambda function with two arguments: file path and Bio.SeqRecord
     - organism    : organism name or NCBI taxid
     - strain      : add this string to the organism name but not to its taxonomy
                     Could be a lambda function with two arguments: file path and Bio.SeqRecord
     - topology    : sequence topology, circular or linear
     - keywords    : keywords to add to annotation
     - ete3        : ete3 object to use when fetching taxonomy
    """

    # Libraries
    import os
    import datetime
    from BCBio import GFF
    import pandas as pd
    from io import StringIO
    from tempfile import NamedTemporaryFile as TempFile
    from collections import OrderedDict
    date = datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime("%d-%b-%Y").upper()

    # Loading taxonomy
    if organism:
        from rotifer.db.ncbi import ncbi
        ncbiobj = ncbi(organism)
        organism = ncbiobj.read("taxonomy", ete3=ete3)

    # Parsing GFF
    h = open(path)
    l = h.readlines()
    h.close()
    n = [ x for x in range(0,len(l)) if l[x] == "##FASTA\n" ][0]
    sio = SeqIO.parse(StringIO("".join(l[n+1:])),"fasta")
    seqrecords = OrderedDict({ x.id: x for x in sio })
    with TempFile() as t:
        # Dump tabular section of GFF to temporary file
        t.write("".join(l[0:n]).encode())
        t.flush()

        # Process each Bio.Seqrecord
        for s in GFF.parse( open( t.name ), base_dict=seqrecords ):
            # Change sequence ID
            seqrecordID = s.id
            if isinstance(seqid,types.FunctionType):
                s.id = seqid(s, path)

            # Annotate date
            s.annotations["date"] = date

            # Annotate taxonomy
            if isinstance(organism, pd.DataFrame):
                s.annotations["organism"] = organism.loc[0,'organism']
                if strain:
                    if isinstance(strain,types.FunctionType):
                        strainStr = strain(path,s)
                    else:
                        strainStr = strain
                    s.annotations["organism"] = s.annotations["organism"] + strainStr
                s.annotations["source"]   = s.annotations["organism"]
                s.annotations["taxonomy"] = organism.loc[0,'taxonomy'].split("; ")

            # Add or set keywords
            if keywords:
                if 'keywords' in s.annotations:
                    s.annotations["keywords"].extend(keywords)
                else:
                    s.annotations["keywords"] = keywords

            # Set description
            if description:
                if isinstance(description,types.FunctionType):
                    descriptionStr = description(path,s)
                else:
                    descriptionStr = description
                s.description = descriptionStr

            # Set topology
            if 'topology' not in s.annotations:
                s.annotations['topology'] = topology

            # Update seqrecord
            if 'gff-version' in s.annotations:
                del(s.annotations['gff-version'])
            if 'sequence-region' in s.annotations:
                del(s.annotations['sequence-region'])
            seqrecords[seqrecordID] = s

        # Close temporary file and return array of seqrecords
        t.close()
        seqrecords = list(seqrecords.values())
        return seqrecords
