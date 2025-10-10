import os
import sys
import numpy as np
import pandas as pd
import pyhmmer as ph
import rotifer
import rotifer.devel.beta.sequence as rdbs
from rotifer.db import ncbi
from rotifer.interval import utils as riu
from rotifer.taxonomy import utils as rtu
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ete3

def digitalize_seqobj(seqobj, alignment=None, msa_name='alignment'):

    """
    Convert sequences from a sequence object into digital format for use with pyhmmer.

    This function takes a sequence object (expected to have a `.df` attribute with
    columns "id" and "sequence"), digitizes the sequences using a aminoacid alphabet,
    and optionally returns a pyhmmer DigitalMSA object.

    Args:
        seqobj: Object
        alignment: bool, optional (default: None)
            If True, returns a DigitalMSA object containing all digitized sequences.
            If False or None, returns a list of pyhmmer DigitalSequence objects.
        msa_name: str, optional (default: 'alignment')
            Name to assign to the DigitalMSA object if `alignment=True`.

    Returns:
        pyhmmer.easel.DigitalMSA or list
            - If `alignment=True`: returns a `DigitalMSA` object containing all
              digitized sequences, suitable for building an HMM.
            - Otherwise: returns a list of `DigitalSequence` objects, one per sequence.

    Example:
        # Digitize sequences without creating a DigitalMSA
        digital_seqs = digitalize_seqobj(seqobj)

        # Create a DigitalMSA for HMM building
        msa = digitalize_seqobj(seqobj, alignment=True, msa_name='my_alignment')
    """

    abc = ph.easel.Alphabet.amino()
    digital_sequences = []
    
    for _, row in seqobj.df.iterrows():
        name = str(row["id"]).encode()
        seq = str(row["sequence"])
        text_seq = ph.easel.TextSequence(name=name, sequence=seq)
        digital_sequences.append(text_seq.digitize(abc))
           
    if alignment:
        for dseq in digital_sequences:
            msa = ph.easel.DigitalMSA(abc, sequences=digital_sequences, name = (msa_name.encode()))    
        return msa
    
    else:
        return digital_sequences
        
def make_hmm(seqobj, hmm_name='alignment', save='alignment.hmm'):
      
    '''
    Build a Hidden Markov Model (HMM) from a sequence object using pyhmmer.
    This function digitizes sequences from a aligned sequence object builds an 
    HMM using pyhmmer's Plan7 Builder, and optionally saves the HMM to disk.

    Args:
        seqobj: Object
        hmm_name: str, optional (default: 'alignment')
            Name assigned to the MSA/HMM. Used internally in the digitalization step.
        save: str or None, optional (default: 'alignment.hmm')
            File path to save the resulting HMM. If None or False, the HMM is not saved.

    Returns:
        tuple of pyhmmer.plan7.HMM
            The HMM object(s) generated from the MSA.

    Example:
        # Build and save an HMM from a sequence object
        hmm = make_hmm(seqobj, hmm_name='my_model', save='my_model.hmm')

        # Build an HMM without saving to disk
        hmm = make_hmm(seqobj, save=None)
    '''
    
    abc = ph.easel.Alphabet.amino()
    msa = digitalize_seqobj(seqobj, alignment=True, msa_name=hmm_name)
    builder = ph.plan7.Builder(abc)
    background = ph.plan7.Background(abc)
    hmm = builder.build_msa(msa, background)
     
    if save:
        with open (save, 'wb') as x:
            hmm[0].write(x)
        print(f'HMM saved in {save}')
        
    return hmm

def hmmscan(sequences, file=None, pfam_database_path='/databases/pfam/Pfam-A.hmm', cpus=0, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to'], rename=True):
    
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
       if hmm_file.is_pressed:
           hmms = list(hmm_file.optimized_profiles())
       else:
           hmms = list(hmm_file)
	
    #Sequences load
    abc = ph.easel.Alphabet.amino()
    if file:
       seqs = ph.easel.SequenceFile(file, digital = True, alphabet=abc)
    
    else:
       if type(sequences) == list:
       		seqobj = rdbs.sequence(sequences)
       		seqs = digitalize_seqobj(seqobj)

       elif type(sequences) == rotifer.devel.beta.sequence.sequence:
       		seqs = digitalize_seqobj(sequences)

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
    
    if rename:
    	df.rename({'aln_target_name': 'sequence',
                   'aln_hmm_name': 'model',
                   'i_evalue': 'evalue',
                   'env_from': 'estart',
                   'env_to': 'eend'}, axis=1, inplace=True)

    return df

def add_arch_to_df(df, column='pid', cpus=0, file=None, pfam_database_path='/databases/pfam/Pfam-A.hmm', inplace=False):
  
    '''
    Add a column pfam with the domain architecture for the input accessions.
    '''
    
    if inplace == False:
    	df = df.copy()
    
    h = hmmscan(df[column].dropna().tolist(), cpus=cpus, file=file, pfam_database_path=pfam_database_path)
    h.rename({'aln_target_name':'sequence','aln_hmm_name':'model','i_evalue':'evalue','env_from':'estart', 'env_to':'eend'}, axis=1, inplace=True)
    arch = riu.filter_nonoverlapping_regions(h.loc[h.groupby(['sequence','model']).score.idxmax()], **riu.config['hmmer']).groupby('sequence').agg(pfam = ('model',lambda x: '+'.join(x.astype(str)))).reset_index()
    arch.rename({'sequence':column}, axis = 1, inplace = True)
    arch = arch.set_index(column).pfam.to_dict()
    df['pfam'] = df[column].map(arch)
    return df
