import os
import sys
import numpy as np
import pandas as pd
import pyhmmer as ph
import rotifer.devel.beta.sequence as rdbs
from rotifer.devel.alpha import trsantos as rdat
from rotifer.db import ncbi
from rotifer.taxonomy import utils as rtu
from rotifer.interval import utils as riu
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ete3

def get_matrix(df, filter_list, rows, columns, n=10, filter_by='pid'):
        filtered_df = df[df[filter_by].isin(filter_list)]
        frequent_rows = filtered_df[rows].value_counts(rows).nlargest(n).index.tolist()
        filtered_df = filtered_df[filtered_df[rows].isin(frequent_rows)]
        matrix = filtered_df.pivot_table(index=rows, columns=columns, aggfunc='size', fill_value=0)
        return matrix

def update_lineage(ndf, preferred="/home/leep/epsoares/projects/databases/data/preferred_taxa.txt"):
    '''
    A function to update the lineage of Gene Neighborhood cursor;
    '''
    tc = ncbi.TaxonomyCursor()
    tax = tc.fetchall(ndf.taxid.dropna().drop_duplicates().tolist())
    tax['taxid'] = tax.taxid.astype(int)
    ndf.classification = ndf.taxid.map(tax.set_index('taxid').classification.to_dict())
    ndf.lineage = rtu.lineage(ndf.classification, preferred_taxa=[ x.strip() for x in open(preferred,"rt")]).tolist()
    return ndf

def extract_envelope(df, seqobj=None, start='estart', end='eend', expand=10, local_database_path='/databases/fadb/nr/nr'):
    '''
    Function to extract the envelope of a model match.
    '''
    if len(df) == 0:
        print(f'Empty input')
        seqobj = rdbs.sequence()

    if seqobj == None:
        seqobj = rdbs.sequence(df.sequence.drop_duplicates().to_list(), local_database_path=local_database_path)

    seqobj.df = seqobj.df.merge(df.rename({'sequence':'id'},axis=1), how='left')
    seqobj.df['end'] = seqobj.df[end].fillna(seqobj.df.length).astype(int)
    seqobj.df['start'] = seqobj.df[start].fillna(1).astype(int)

    if expand:
        seqobj.df['start'] = np.where(seqobj.df.start<=expand,1,seqobj.df.start-expand+1)
        seqobj.df['end'] = np.where(seqobj.df.end+expand > seqobj.df.length, seqobj.df.length, seqobj.df.end+expand)

    seqobj.df.sequence = seqobj.df.apply(lambda row: row['sequence'][row['start'] - 1:row['end']], axis=1)
    seqobj.df['pid'] = seqobj.df['id']
    seqobj.df['id'] = seqobj.df['id'] + '/' + seqobj.df[start].astype(str) + '-' + seqobj.df[end].astype(str)
    seqobj._reset()

    return seqobj

def to_network(df, target=['pfam'], ftype=['CDS'], interaction=True, ignore = [], strand = True):
    if isinstance(ftype, str):
        ftype = [ftype]
    if isinstance(target, str):
        itarget = [target]

    w = df.copy()
    if strand:
        w = df.neighbors(strand='same')
        w['rid'] = list(range(1,len(w)+1))
        w.rid = w.rid * w.strand
        w.sort_values(['rid'], inplace=True)

    # Building the source column
    w = w.query(f'type == "{ftype}"').block_id.reset_index().drop('index', axis=1)
    w['source'] = df[target[0]]
    for col in target[1:]:
        w['source'] = np.where(w['source'].isna(), df[df.type == ftype][col], w['source'])
    w['source'] = w['source'].fillna("?").str.split('+')
    w = w.explode(column='source')
    if ignore:
        w = w[~w.source.isin(ignore)].copy()

    # Building target data
    w['tblock_id'] = w['block_id'].shift(-1)
    w['target'] = w['source'].shift(-1)

    # Selecting same block rows
    w = w[w.block_id == w.tblock_id.shift(1)].copy()

    # Fix source target order when not restricted to the same strand
    sameprotein = (w.index.to_series() == w.index.to_series().shift(1))
    if not strand:
        reverse = (~sameprotein) & (w.source > w.target)
        w.loc[reverse,['source','target']] = w.loc[reverse, ['target', 'source']].values

    if interaction:
        w['interaction'] = np.where(sameprotein, 'fusion', 'neighbor')
        w = w.groupby(['source', 'target', 'interaction'])
    else:
        w = w.groupby(['source', 'target'])

    w = w.agg(weight=('block_id', 'count'), blocks=('block_id', 'nunique')).reset_index()
    return w

def make_palette(categories):
    """
    Generate a color-blind safe palette dictionary 
    for the given list of category names (max 10 entries).
    """
    safe_colors = [
        "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#999999",
        "#117733", "#882255"
    ]

    # Keep only as many as the safe palette allows
    categories = categories[:len(safe_colors)]

    palette = {cat: safe_colors[i] for i, cat in enumerate(categories)}

    return palette

def attribute_table(
    ndf,
    columns_list=['pid', 'pfam', 'compact_total', 'compact_same_strand', 'kingdom', 'phylum', 'class', 'classification'],
    tax_parser= True,
    color_by_tax='phylum',
    number_of_taxa=5,
    save=None):

    '''
    Make a attribute table using the Neighborhood dataframe with a architecure column, 
    select and reorder columns to match those required by TreeViewer (leaf identifier 
    as first column). Also permits automatic coloring by taxonomy.
    '''

    ndfc = ndf.compact()
    ndfcs = ndf.select_neighbors(strand = True).compact()
    att = ndf[ndf.pid.isin(ndf.query('query == 1').pid)].drop_duplicates('pid')
    att['compact_total'] = att.block_id.map(ndfc['compact'].to_dict())
    att['compact_same_strand'] = att.block_id.map(ndfcs['compact'].to_dict())

    if tax_parser:
        tax = rdat.taxon_summary(att).set_index('taxid')
        for col in tax.columns:
            att[col] = att['taxid'].map(tax[col].to_dict())

    if columns_list:
        att = att[columns_list]

    if color_by_tax:
        palette = make_palette(att[color_by_tax].value_counts().head(number_of_taxa).index.tolist())
        att['Color'] = att[color_by_tax].map(palette).fillna('#000000')

    if save:
        att.to_csv(save, sep = '\t', index = False)
        print(f'Table saved to {save}')

    else:
        return att

def make_heatmap(
        ndf,
        name = None,
        domain_list = [],
        format_table = True,
        color_list = ['white','blue','green','red'],
        cbar = False,
        fmt = 'd',
        annot = False,
        linewidths = 1,
        tree_file = None,
        tree = False,
        figsize=None
    ):
    '''
    Creates a heatmap from a dataframe output of full annotate
    to represent main occurrence. The colors are correspondent 
    to the values in the dataframe, by exemple, four colors 
    corresponds to the values: 0, 1, 2 and 3, respectively. 
    The columns in the dataframe will be the columns in the 
    heatmap, and each line represents a occurrence from the 
    respective domain and the colors the type of correspondence, 
    that is determined by the numbers utilized. 0 equals to non 
    correlation, 1 correlation by neighborhood, 2 correlation by 
    fusions and 3 by both. If a tree file is delivered by the user
    the heatmap will sort the occurrences by the order of the tree.
    '''
 
    if format_table == True:
        hm = ndf.set_index('pid')[['block_id']]
        hm['neighbors'] = hm.block_id.map(ndf.query('query == 0').groupby('block_id').agg(neighbors = ('pfam','sum')).neighbors.to_dict())
        hm['pfam'] = hm.block_id.map(ndf.query('query == 1').set_index('block_id').pfam.to_dict())
        for x in domain_list:
            inpfam = hm['pfam'].str.contains(x, na = False)
            inneighbors = hm['neighbors'].str.contains(x, na = False)
            hm[f'{x}'] = ((inpfam & inneighbors).astype(int)*3 + (inpfam & ~inneighbors).astype(int)*2 + (~inpfam & inneighbors).astype(int)*1 + (~inpfam & ~inneighbors).astype(int)*0)
        df = hm[domain_list]
                                                                   
    if tree == True:
        t = ete3.Tree(f'{tree_file}')
        df = df.reindex(df.reindex(t.get_leaf_names()).index.str.replace('\'','')).fillna(0).astype(int)
    
    if figsize:
        figsize=figsize
    else:
        figsize=(len(df.columns),len(df)/5)
  
    cmap = colors.ListedColormap(color_list)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(data = df, cmap = cmap, cbar = cbar, fmt = fmt, annot = annot, linewidths = linewidths, ax=ax)
    plt.savefig(f'{name}', bbox_inches='tight')
    plt.close(fig)

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

def add_arch_to_df(df, column='pid', cpus=0, file=None, pfam_database_path='/databases/pfam/Pfam-A.hmm'):
  
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
      
    
  
