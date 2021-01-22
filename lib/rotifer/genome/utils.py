#!/usr/bin/env python3

# Core libraries
import sys
import types

# Pandas
import numpy as np
import pandas as pd

# Rotifer libraries

def seqrecord2df(seqrecs, exclude_type=[], autopid=False, assembly=None, codontable='Bacterial', block_id=-1):
    '''
    Extract BioPython's SeqRecord features data to a Pandas dataframe
    Arguments:
      seqrecs      : rotifer.genome.io.parse generator or Bio.SeqRecord object(s)
      exclude_type : exclude features by type

      autopid      : auto-generate missing PIDs from locus_tag
       You can set autopid to:
        - 'auto', for transforming the locus_tag into unique PIDs, or
        - 'copy', to copy the locus_tag, if there is no alternative splicing

      codontable   : Genetic code name for translating CDSs
    '''
    import Bio
    from math import ceil, log10
    from Bio.Alphabet import generic_dna
    from rotifer.genome.data import NeighborhoodDF

    # Make sure first argument is a list
    dfstack = []
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
                    try:
                        plen = len(ft.translate(seqrecord, table=codontable))
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
            gene = qualifiers['gene'][0] if 'gene' in qualifiers else '.'
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
            if len(ft.location.parts) == 1:
                start = ft.location.start+1
                end = int(ft.location.end)
                location = [[start, end]]
            else:
                location = []
                for l in ft.location.parts:
                    start = l.start+1
                    end = int(l.end)
                    if location:
                        if ((end - location[-1][1]) * strand) <= 0:
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
                    'block_id': block_id,
                    'rid':0,
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
                    'classification': taxonomy,
                    'feature_order':feature_order[feature_type],
                    'internal_id':internal_id})

            # Increment feature counters
            feature_order[feature_type] += 1
            internal_id += 1

        # Decrement block_id
        block_id -= 1

    # Build dataframe and return
    if len(data) > 0:
        df = pd.DataFrame(data)
        if exclude_type:
            df = df[~df['type'].isin(exclude_type)]
        df = NeighborhoodDF(df, update_lineage='classification')
        return(df)
    else:
        return(pd.DataFrame())

