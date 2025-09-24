#!/usr/bin/env python3

# Core libraries
import sys
import types

# Pandas
import pandas as pd
import numpy as np

# Rotifer libraries
import rotifer
logger = rotifer.logging.getLogger(__name__)

# Data
_columns = ['nucleotide', 'start', 'end', 'strand', 'nlen', 'block_id', 'query', 'pid', 'type', 'plen', 'locus', 'seq_type', 'assembly', 'gene', 'origin', 'topology', 'product', 'taxid', 'organism', 'lineage', 'classification', 'feature_order', 'internal_id', 'is_fragment']

def seqrecords_to_dataframe(seqrecs=None, exclude_type=[], autopid=False, assembly=None, codontable='Bacterial'):
    '''
    Extract BioPython's SeqRecord features data to a Pandas dataframe
    Arguments:
      seqrecs      : rotifer.genome.io.parse generator, Bio.SeqIO generator
                     or (list of) Bio.SeqRecord object(s)
      exclude_type : exclude features by type
      autopid      : auto-generate missing PIDs from locus_tag

          You can set autopid to:
          - 'auto', for transforming the locus_tag into unique PIDs, or
          - 'copy', copy the locus_tag if there is no alternative splicing

      codontable   : Genetic code name for translating CDSs
    '''
    import re
    import os
    import Bio
    from math import ceil, log10
    from rotifer.genome.data import NeighborhoodDF
    ncbiaccre = re.compile(r'^GC[AF]_[0-9]+\.[0-9]+') # Regular expression to extract NCBI assembly IDs

    # Make sure first argument is a list
    if not (isinstance(seqrecs,types.GeneratorType) or isinstance(seqrecs,list) or isinstance(seqrecs,Bio.SeqIO.Interfaces.SequenceIterator)):
        seqrecs = [ seqrecs ]

    # Process each seqrecord
    data = []
    for seqrecord in seqrecs:
        # Check input
        if not isinstance(seqrecord, Bio.SeqRecord.SeqRecord):
            logger.error(f'Element is not a seqrecord: {seqrecord}')
            continue
        if not seqrecord.features:
            #logger.error(f'Bio.Reqrecord {seqrecord.id} contains has features to extract')
            continue
        digits = max(6,ceil(log10(len(seqrecord.features))))
        assemblyID = None

        # Extract SeqRecord data
        nlen = len(seqrecord)
        if seqrecord.annotations:
            annotations = seqrecord.annotations
            topology = annotations['topology'] if 'topology' in annotations else 'linear'
            organism = annotations['organism'] if 'organism' in annotations else 'Unknown'
            taxonomy = '; '.join(annotations['taxonomy']) if 'taxonomy' in annotations else ''
            for x in seqrecord.dbxrefs:
                if 'Assembly:' in x:
                    assemblyID = x.split(':')[-1]
        else:
            topology = 'linear'
            organism = 'Unknown'
            taxonomy = pd.NA

        # SeqRecord doesn't contain a reference to its assembly ID
        # Let's assume all sequences in the current file belong to
        # the same assembly
        if assemblyID == None:
            if hasattr(seqrecs,"stream"): # Bio.SeqIO.Interfaces.SequenceIterator!
                if hasattr(seqrecs.stream,"filename") and callable(seqrecs.stream.filename):
                    assemblyID = seqrecs.stream.filename() # rotifer.io.fileinput.FileInput
                elif hasattr(seqrecs.stream,"name"):
                    assemblyID = seqrecs.stream.name # Another _io class
                if assemblyID != None:
                    try:
                        assemblyID = os.path.basename(assemblyID)
                        if ncbiaccre.match(assemblyID):
                            assemblyID = ncbiaccre.search(assemblyID).group(0)
                    except:
                        logger.debug(f"Assembly {str(assemblyID)} doesn't parse as a string")
        if assemblyID == None:
            if assembly:
                assemblyID = assembly
            else:
                assemblyID = seqrecord.id

        # Process source feature
        taxid = np.nan
        seq_type = 'chromosome'
        ftsource = [ x for x in seqrecord.features if x.type == "source" ]
        if len(ftsource):
            ftsource = ftsource[0]
            if "db_xref" in ftsource.qualifiers:
                for dbx in ftsource.qualifiers['db_xref']:
                    if "taxon:" == dbx[0:6]:
                        taxid = int(dbx[6:])
            if 'plasmid' in ftsource.qualifiers:
                seq_type = 'plasmid'

        # Extract data for each feature
        feature_order = {}
        internal_id = -1
        seqrecord.features.sort(key=lambda x: (x.location.start, x.location.end))
        for ft in seqrecord.features:
            internal_id += 1 # Internal ID depends solely on the features annotated in the data source
            qualifiers = ft.qualifiers

            # Feature type
            feature_type = ft.type

            # Pseudogenes
            if feature_type not in feature_order:
                feature_order[feature_type] = 0
            if feature_type == 'pseudogene':
                feature_type = 'PSE'
                if 'PSE' not in feature_order:
                    feature_order['PSE'] = 0
            if exclude_type and feature_type in exclude_type:
                continue

            # PID
            locus = qualifiers['locus_tag'][0] if 'locus_tag' in qualifiers else seqrecord.id + f'.{internal_id:0{digits}}'
            plen = sum([ len(x) for x in ft.location.parts ])
            pid = np.nan
            if feature_type == 'CDS':
                if 'pseudo' in qualifiers:
                    feature_type = 'PSE'
                    if 'PSE' not in feature_order:
                        feature_order['PSE'] = 0
                else:
                    if 'translation' in qualifiers:
                        plen = len(qualifiers['translation'][0])
                    else:
                        selectedTable = codontable
                        if 'transl_table' in qualifiers:
                            selectedTable = int(qualifiers['transl_table'][0])
                        try:
                            plen = len(ft.translate(seqrecord, table=selectedTable, cds=False, to_stop=True))
                        except:
                            feature_type = 'PSE'
                            if 'PSE' not in feature_order:
                                feature_order['PSE'] = 0
                if feature_type == 'CDS':
                    if 'protein_id' in qualifiers:
                        pid = qualifiers['protein_id'][0]
                    elif 'ID' in qualifiers:
                        pid = qualifiers['ID'][0]
                    elif autopid == 'auto':
                        pid = f'{locus}.p{feature_order["CDS"]:0{digits}}'
                    elif autopid == 'copy':
                        pid = locus

            # Other feature attributes
            gene = qualifiers['gene'][0] if 'gene' in qualifiers else np.nan
            product = np.nan
            for tag in ['product','inference']:
                if tag in qualifiers:
                    product = qualifiers[tag][0]
                    break

            # Strand
            try:
                strand = int(ft.location.strand)
            except:
                strand = 1

            # Feature location: check if multiple locations imply running through the origin
            origin = 0
            location = []
            for l in ft.location.parts:
                start = l.start+1
                end = int(l.end)
                if location:
                    if ((end - location[-1][1]) * strand) <= 0:
                        if end == nlen and location[-1][0] == 1:
                            location.insert(-1,[start,end])
                        else:
                            location.append([start,end])
                        origin = 1
                    else:
                        location[-1][0] = min(location[-1][0],start)
                        location[-1][1] = max(location[-1][1],end)
                else:
                    location.append([start,end])

            # Store each location
            for l in location:
                data.append({
                    'nucleotide': seqrecord.id,
                    'start': l[0],
                    'end': l[1],
                    'strand': strand,
                    'nlen': nlen,
                    'block_id': seqrecord.id,
                    'query':0,
                    'pid': pid,
                    'type': feature_type,
                    'plen': plen,
                    'locus': locus,
                    'seq_type': seq_type,
                    'assembly': assemblyID,
                    'gene':gene,
                    'origin':origin,
                    'topology': topology,
                    'product': product,
                    'taxid': taxid,
                    'organism': organism,
                    'lineage':pd.NA,
                    'classification': taxonomy,
                    'feature_order':feature_order[feature_type],
                    'internal_id':internal_id,
                    'is_fragment':0})

            # Increment feature counters
            feature_order[feature_type] += 1

        #END: for ft in seqrecord.features
    #END: for seqrecord in seqrecs

    # Build dataframe and return
    if len(data) > 0:
        df = pd.DataFrame(data)
        df = NeighborhoodDF(df, update_lineage='classification')
        return(df)
    else:
        return(NeighborhoodDF(pd.DataFrame(columns=_columns)))

def gembase(ndf, name, strain=0, inplace=True, protein_dict=None):
    from datetime import datetime
    isProtein = ndf.type == "CDS"
    assemblyChange = ndf.assembly != ndf.assembly.shift(1)
    gembase_strain = (assemblyChange.cumsum() + int(strain) - 1).astype(str).str.zfill(5)
    protein_number = isProtein.cumsum() - pd.Series(np.where(assemblyChange, isProtein.cumsum() - 1, np.NaN)).ffill()
    protein_number = protein_number.astype(int).astype(str).str.zfill(5)
    protein_number = pd.Series(np.where(isProtein, protein_number, ""))
    if protein_dict:
        gembase = gembase_strain + ":" + protein_number
        gembase = gembase.map(protein_dict).fillna("")
    else:
        contigChange = ndf.nucleotide != ndf.nucleotide.shift(1)
        ncontig = contigChange.cumsum()
        contig_number = np.where(assemblyChange, ncontig - 1, np.NaN)
        contig_number = ncontig - pd.Series(contig_number).ffill()
        contig_number = contig_number.astype(int).astype(str).str.zfill(5)
        contig_place = pd.Series(np.where(contigChange | (ndf.nucleotide != ndf.nucleotide.shift(-1)), "b", "i"))
        gembase = name + "." + datetime.now().strftime("%m") + datetime.now().strftime("%y")
        gembase = gembase + "." + gembase_strain
        gembase = pd.Series(np.where(isProtein, gembase + "." + contig_number + contig_place + "_", ""))
    if inplace == True:
        ndf['gembase'] = gembase
    else:
        return gembase

