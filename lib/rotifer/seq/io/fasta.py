#!/usr/bin/env python3

from io import StringIO as iStringIO
from Bio import SeqIO
import pandas as pd
from rotifer.seq.alignment import _seq2matrix
class fasta:
    def __init__(self, alignment,input_format, **kwargs):
        self.alignment = alignment
        self.input_format = input_format

    def run(self):
        idd=[]
        seq=[]
        fasta_sequences = SeqIO.parse(alignment, input_format)

        for fasta in fasta_sequences:
            idd.append(fasta.id)
            seq.append(str(fasta.seq))

        data={'sequence_id':idd, 'sequence':seq}

        df = pd.DataFrame(data=data)
        df['sequence'] = df['sequence'].str.ljust(df['sequence'].str.len().max(), '-')
        metadata = {'sequence_id': 'sequence_id',
                    'sequence': 'sequence'}

        ma = _seq2matrix(df[['sequence']], metadata = metadata)
        df = df.drop(columns = ['sequence'])

        tuple_columns = ()

        for col in df.columns:
            if col == 'sequence_id':
                tuple_columns += ('_sequence_id', col)
            else:
                tuple_columns += ('other_information', col)

        df.columns = pd.MultiIndex.from_tuples(
                                        [tuple_columns]
                                        )

        df = pd.concat([df,ma], axis = 1)
        return df, metadata


        #return cls(df,columns = df.columns, metadata = metadata)
