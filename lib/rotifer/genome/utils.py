#!/usr/bin/env python3

# Core libraries
import sys
import types

# Pandas
import pandas as pd

# Rotifer libraries

# Data
_columns = ['nucleotide', 'start', 'end', 'strand', 'nlen', 'block_id', 'query', 'pid', 'type', 'plen', 'locus', 'seq_type', 'assembly', 'gene', 'origin', 'topology', 'product', 'organism', 'lineage', 'classification', 'feature_order', 'internal_id', 'is_fragment']

def seqrecords_to_dataframe(seqrecs, exclude_type=[], autopid=False, assembly=None, codontable='Bacterial', block_id=-1):
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
      block_id     : initial value of the default block_id
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
            print(f'Element is not a seqrecord: {seqrecord}', file=sys.stderr)
            continue
        if not seqrecord.features:
            #print(f'Bio.Reqrecord {seqrecord.id} contains has features to extract', file=sys.stderr)
            continue
        digits = max(6,ceil(log10(len(seqrecord.features))))

        # Extract SeqRecord data
        nlen = len(seqrecord)
        if seqrecord.annotations:
            annotations = seqrecord.annotations
            topology = annotations['topology'] if 'topology' in annotations else 'linear'
            organism = annotations['organism'] if 'organism' in annotations else 'Unknown'
            taxonomy = '; '.join(annotations['taxonomy']) if 'taxonomy' in annotations else ''
            assemblyID = [ x.split(':')[-1] for x in seqrecord.dbxrefs if 'Assembly:' in x ]
            if assemblyID:
                assembly = assemblyID[0]
        else:
            topology = 'linear'
            organism = 'Unknown'
            taxonomy = pd.NA

        # SeqRecord doesn't contain a reference to its assembly ID
        # Let's assume all sequences in the current file belong to
        # the same assembly
        if assembly == None:
            if hasattr(seqrecs,"stream"): # Bio.SeqIO.Interfaces.SequenceIterator!
                if hasattr(seqrecs.stream,"filename") and callable(seqrecs.stream.filename):
                    assembly = seqrecs.stream.filename() # rotifer.io.fileinput.FileInput
                elif hasattr(seqrecs.stream,"name"):
                    assembly = seqrecs.stream.name # Another _io class
                if assembly != None:
                    assembly = os.path.basename(assembly)
                    if ncbiaccre.match(assembly):
                        assembly = ncbiaccre.search(assembly).group(0)
        if assembly == None:
            assembly = seqrecord.id

        # Extract data for each feature
        seq_type = 'chromosome'
        feature_order = {}
        internal_id = 0
        seqrecord.features.sort(key=lambda x: (x.location.start, x.location.end))
        for ft in seqrecord.features:
            qualifiers = ft.qualifiers

            # Feature type
            feature_type = ft.type
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
            pid = ''
            if feature_type == 'CDS':
                if 'pseudo' in qualifiers:
                    feature_type = 'PSE'
                    if 'PSE' not in feature_order:
                        feature_order['PSE'] = 0
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
                    elif autopid == 'auto':
                        pid = f'{locus}.p{feature_order["CDS"]:0{digits}}'
                    elif autopid == 'copy':
                        pid = locus

            # Other feature attributes
            if 'plasmid' in qualifiers:
                seq_type = 'plasmid'
            gene = qualifiers['gene'][0] if 'gene' in qualifiers else ''
            product = ''
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
                    'block_id': block_id,
                    'query':0,
                    'pid': pid,
                    'type': feature_type,
                    'plen': plen,
                    'locus': locus,
                    'seq_type': seq_type,
                    'assembly': assembly,
                    'gene':gene,
                    'origin':origin,
                    'topology': topology,
                    'product': product,
                    'organism': organism,
                    'lineage':pd.NA,
                    'classification': taxonomy,
                    'feature_order':feature_order[feature_type],
                    'internal_id':internal_id,
                    'is_fragment':0})
                internal_id += 1

            # Increment feature counters
            feature_order[feature_type] += 1

        # Decrement block_id and last feature id
        block_id -= 1
        #END: for ft in seqrecord.features
    #END: for seqrecord in seqrecs

    # Build dataframe and return
    if len(data) > 0:
        df = pd.DataFrame(data)
        df = NeighborhoodDF(df, update_lineage='classification')
        return(df)
    else:
        return(NeighborhoodDF(pd.DataFrame(columns=_columns)))

