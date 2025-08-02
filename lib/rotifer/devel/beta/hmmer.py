import os
import sys
import numpy as np
import pandas as pd
import pyhmmer as ph
import rotifer.devel.beta.sequence as rdbs
from rotifer.db import ncbi
from rotifer.interval import utils as riu
from rotifer.taxonomy import utils as rtu
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ete3

def hmmscan(sequences, file=None, pfam_database_path='/databases/pfam/Pfam-A.hmm', cpus=0, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to']):
    
    '''
    Perform an hmmscan of protein sequences against a Pfam HMM database.
    This function accepts an sequence object or a FASTA file path, and 
    returns a pandas DataFrame summarizing domain hits and alignment attributes.
    ----------
    Parameters
    ----------
    sequences: sequence object
    file: str, optional
        Path to a FASTA file containing the query sequences. If provided, `sequences`
        is ignored and sequences are read from this file directly.
    pfam_database_path : str, optional
        Filesystem path to the Pfam-A HMM database file.
    cpus: int, optional
        Number of CPU threads to allocate for the hmmscan search. A value of 0
        lets HMMER autodetect available cores.
    columns: list of str, optional
        Subset of result column names to include in the output DataFrame.
        Default columns include basic domain and alignment metrics.
    '''

    #HMM load
    with ph.plan7.HMMFile(pfam_database_path) as hmm_file:
        hmms = list(hmm_file.optimized_profiles())

    #Sequences load
    abc = ph.easel.Alphabet.amino()
    if file:
        seqs = ph.easel.SequenceFile(file, digital = True, alphabet=abc)
    else:
        if type(sequences) == list:
            sequences = rdbs.sequence(sequences)
        sequences.to_file('sequences.fasta')
        seqs = ph.easel.SequenceFile('sequences.fasta', digital=True, alphabet=abc)
        os.remove('sequences.fasta')

    #Hmmscan run and file processment
    h = list(ph.hmmer.hmmscan(seqs, hmms, cpus=cpus))
    r = []
    for th in h:
        for x in th:
            for y in x.domains:
                r.append({
                    # pyhmmer.plan7.Domain attributes
                    "hit":                   y.hit,
                    "bias":                  y.bias,
                    "c_evalue":              y.c_evalue,
                    "correction":            y.correction,
                    "env_from":              y.env_from,
                    "env_to":                y.env_to,
                    "env_score":             y.envelope_score,
                    "i_evalue":              y.i_evalue,
                    "pvalue":                y.pvalue,
                    "score":                 y.score,

                    # pyhmmer.plan7.Alignment attributes
                    "aln_domain":            y.alignment.domain,
                    "aln_hmm_accession":     y.alignment.hmm_accession.decode(),
                    "aln_hmm_from":          y.alignment.hmm_from,
                    "aln_hmm_name":          y.alignment.hmm_name.decode(),
                    "aln_hmm_sequence":      y.alignment.hmm_sequence,
                    "aln_hmm_to":            y.alignment.hmm_to,
                    "aln_hmm_length":        y.alignment.hmm_length,
                    "aln_identity_sequence": y.alignment.identity_sequence,
                    "aln_target_from":       y.alignment.target_from,
                    "aln_target_name":       y.alignment.target_name.decode(),
                    "aln_target_sequence":   y.alignment.target_sequence,
                    "aln_target_to":         y.alignment.target_to,
                    'aln_target_length':     y.alignment.target_length
                    })
    df = pd.DataFrame(r)

    if columns:
        df = df[columns]

    return df

def add_arch_to_df(df, column='pid'):
    '''
    Add a column pfam with the domain architecture for the input accessions.
    '''
    h = hmmscan(df[column].dropna().tolist())
    h.rename({'aln_target_name':'sequence','aln_hmm_name':'model','i_evalue':'evalue','env_from':'estart', 'env_to':'eend'}, axis=1, inplace=True)
    arch = riu.filter_nonoverlapping_regions(h, **riu.config['hmmer']).groupby('sequence').agg(pfam = ('model',lambda x: '+'.join(x.astype(str)))).reset_index()
    arch.rename({'sequence':column}, axis = 1, inplace = True)
    arch = arch.set_index(column).pfam.to_dict()
    df['pfam'] = df[column].map(arch)
    return df
