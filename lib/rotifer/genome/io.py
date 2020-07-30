#!/usr/bin/env python3

# Core libraries
import types
import itertools

# Pandas
import numpy as np
import pandas as pd

# Biopython
from Bio import SeqIO
from Bio.Alphabet import generic_dna

# Rotifer libraries
from rotifer.core.functions import not_kwargs
from rotifer.taxonomy.utils import lineage as rtlineage

def parse(handle, informat, *args, **kwargs):
    """
    Load genomic data from sequence file(s).
    This method crates an iterator capable of return SeqRecords.

    Arguments:
     - handle   : one or more file handles or filenames
     - informat : either gff or a Bio.SeqIO supported format
     - retvalue : object(s) type to return. Supported objects are

     Other arguments depend on the parser (see methods below)
    """
    if not isinstance(handle,list):
        handle = [handle]
    for eachhandle in handle:
        if informat == "gff":
            for s in gff(eachhandle, *args, **kwargs):
                yield s
        else:
            for s in SeqIO.parse(eachhandle,informat,alphabet=generic_dna):
                yield s

def gff(path, seqid=None, assembly=None, organism=None, strain=None, description=None, topology='linear', keywords=[], ete3=None):
    """
    Parse GFF file

    Arguments:
     - path        : file path or file handle
     - seqid       : rename sequence identifiers using a user-defined function
                     This function will be given two arguments:
                      - the sequence object (a Bio.SeqRecord.SeqRecord object)
                      - the input file path
     - description : set sequence description
                     Could be a lambda function to apply to the file path
     - assembly    : genome assembly accession
                     Could be a lambda function to apply to the file path
     - organism    : organism name or NCBI taxid
     - strain      : add this string to the organism name but not to its taxonomy
                     Could be a lambda function to apply to the file path
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
    sio = SeqIO.parse(StringIO("".join(l[n+1:])),"fasta", alphabet=generic_dna)
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
                        strainStr = strain(path)
                    else:
                        strainStr = strain
                    s.annotations["organism"] = s.annotations["organism"] + strainStr
                s.annotations["source"]   = s.annotations["organism"]
                s.annotations["taxonomy"] = organism.loc[0,'taxonomy'].split("; ")

            # Fix assembly name
            if assembly:
                if isinstance(assembly,types.FunctionType):
                    assemblyStr = assembly(path)
                else:
                    assemblyStr = assembly
                s.dbxrefs.append(f'Assembly:{assemblyStr}')
            else:
                s.dbxrefs.append(f'Assembly:{s.name}')

            # Add or set keywords
            if keywords:
                if 'keywords' in s.annotations:
                    s.annotations["keywords"].extend(keywords)
                else:
                    s.annotations["keywords"] = keywords

            # Set description
            if description:
                if isinstance(description,types.FunctionType):
                    descriptionStr = description(path)
                else:
                    descriptionStr = assembly
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
