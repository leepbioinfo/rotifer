def load_seq_scan(name, folder):
    '''
    load a seqscan result into a dataframe
    '''
    pre = f'{folder}/{name}'
    info = pd.read_csv(f'{pre}.c100i100.tsv', sep='\t', names=['c100i100', 'pid'])
    info = info.merge(pd.read_csv(f'{pre}.c80i70.tsv', sep='\t', names=['c80i70', 'c100i100']), how="left")
    info = info.merge(pd.read_csv(f'{pre}.c80e3.tsv', sep='\t', names=['c80e3', 'c80i70']), how="left")
    info = info.merge(pd.read_csv(f'{pre}.aravind.scan.arch', sep='\t', names=['c100i100', 'aravind'], usecols=[0,1], skiprows=[0]), how="left")
    info = info.merge(pd.read_csv(f'{pre}.pfam.hmmscan.arch', sep='\t', names=['c100i100', 'pfam'], usecols=[0,1], skiprows=[0]), how="left")
    return info


def cluster2aln(group_cluster,df,esl_index_file, grouper='c80e3', redundancy_cluster='c80i70', fast=True, query=False):
    import tempfile
    from subprocess import Popen,PIPE
    import os

    with tempfile.TemporaryDirectory() as tmpdirname:
        if query:
            df.query(query)[redundancy_cluster].drop_duplicates().dropna().to_csv(f'{tmpdirname}/accs', index=None, header=None)
        else:
            df[df[grouper]  == group_cluster][redundancy_cluster].drop_duplicates().dropna().to_csv(f'{tmpdirname}/accs', index=None, header=None)
        if not os.path.exists(esl_index_file+'.ssi'):
            Popen(f'esl-sfetch --index {esl_index_file}',stdout=PIPE, shell=True).communicate()
        Popen(f'esl-sfetch -f {esl_index_file} {tmpdirname}/accs > {tmpdirname}/accs.fa',stdout=PIPE, shell=True).communicate()
        b = sequence(f'{tmpdirname}/accs.fa').realign(fast=fast)
        return b


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def cluster_Co_occurrence(df, count='c80e3', freq_cutoff = 0.3, only_query=True, annotation=False):
    ''' Function to count the co-occurence of clusters within the NeighborhoodDF.
    The count parameter shoul use to define the clster that would be analysed and 
    the freq_cutoff is to use to define a minimum cut_off to display
    '''
    x = df[['block_id', count]].merge(df[['block_id', count]], how='outer', on='block_id')
    x.columns = ['block_id', 'query_cluster', 'neighbor_cluster']
    x = x[x.query_cluster != x.neighbor_cluster]
    if only_query:
        xx = x[x['query_cluster'].isin(df.query('query ==1').pid.unique())]
    else:
        xx = x

    xxx = xx.groupby(['query_cluster', 'neighbor_cluster']).block_id.nunique().reset_index()
    xxxx = xx.groupby('query_cluster').block_id.nunique().rename('query_blocks').reset_index().merge(xxx, how='left').sort_values(['query_blocks', 'block_id'], ascending=False)
    xxxx['query_freq'] = xxxx.block_id/xxxx.query_blocks
    if not annotation:
       return xxxx.query('query_freq >= @freq_cutoff')

    andf= df.groupby(count).agg(pfam = ('pfam', count_series), aravind = ('aravind', count_series)).reset_index()
    xxxx = xxxx.merge(andf.rename({count:'query_cluster', 'pfam':'query_pfam', 'aravind':'query_aravind'}, axis=1), how='left')
    xxxx = xxxx.merge(andf.rename({count:'neighbor_cluster', 'pfam':'neighbor_pfam', 'aravind':'neighbor_aravind'}, axis=1), how='left')
    return xxxx.query('query_freq >= @freq_cutoff')



def count_series(series, normalize=False, cut_off=False, count='domain'):
    ''' 
    Function to flatten a pd.Series.value_counts and count the number of domains recognized in the whole proteins set of proteins
    Mind that by defaul, it count domain not architecture. If one wants to count architecture, it should change the count option to 'architecture'.
    The normalize options is to write the results as frequency.
    The cut_off options is to print results only above a given threashold. If the normalize function is True, you should give the cutoff value as pct. 
    Usage:
    in: count_series(df.pfam, normalize=True, cut_off=10)
    out:'T6SS_HCP(15.15%), SP(11.36%)'
    in: df.groupby('c80e3').agg({'pfam':count_series})
                                                 pfam
    c80e3                                            
    AAN66278.1                                       
    ABO22067.1                           SP(1), TM(1)
    ABU79335.1  Gal_mutarotas_2(1), Glyco_hydro_31(1)
    ABV15410.1                           SP(1), TM(1)
    ABV15414.1                                  TM(1)

    '''
    l = []
    if count == 'architecture':
        s = series.dropna().value_counts(normalize=normalize)
    else:
        s = series.str.split('+').explode().dropna().value_counts(normalize=normalize)

    if cut_off:
        cut_off = cut_off/100
        s = s.where(lambda x: x>= cut_off).dropna()
    for y,z in s.iteritems():
        if normalize:
            l.append(f'{y}({100 * z:.2f}%)')
        else:
            l.append(f'{y}({z})')
    return ', '.join(l)

def fetch_seq(seqs):
    import time
    from Bio import Entrez
    import rotifer.db.ncbi 
    from rotifer.db.ncbi import NcbiConfig
    from rotifer.devel.beta.sequence import sequence
    from rotifer.devel.alpha.gian_func import chunks
    if isinstance(seqs, list):
        seqs = chunks(seqs, 200)
    else:
     return print('add seuquences as list object')

    seq_string = ''
    for x in seqs:
        f = Entrez.efetch(db = 'protein', rettype = 'fasta', retmode = 'text', id = ','.join(x), api_key = NcbiConfig['api_key']).read()
        seq_string = seq_string + f
        time.sleep(1)
    return sequence.from_string(seq_string)

def annotation(seqobj, start, end, annotation):
    ''' I sill have to finish it, but I  am almost done'''
    s = seqobj.copy()
    t = pd.Series(list(s.df.iloc[0,1]))
    t.iloc[:] =' '
    size = end - start
    size_an = len(annotation)
    m = int((start+end)/2)
    t.iloc[start:end] ='∞'
    t.iloc[start] ='∞'
    t.iloc[end] ='∞'
    t.iloc[m] ='B'
    t.iloc[m+1] ='1'
    s.df = pd.concat([pd.DataFrame([['teste', ''.join(t.to_list()), 'annotation']], columns=['id', 'sequence', 'type']),s.df])
    return s.view()

def hmmsearch_full2pandas (file, error_lines = True):
    """
    Function to load the result of hmmsearch or scan ang obtain the results for the full protein
    If error_lines = False, it skip the results that gave error loading to the dataframe.
    Those errors only happens if any field present an string with more than two spaces in sequence (mostly found in  descrption field )
    """
    import re
    from io import StringIO
    
    with open(file, 'r') as f:
        text = f.read()
    columns = 'Evalue  score  bias    BD_Evalue  BD_score  BD_bias    exp  N  Sequence       Description'.split()
    u = '------- ------ -----    ------- ------ -----   ---- --  --------       -----------'
    l = 'Domain annotation for each sequence'
    match = re.findall(f'{u}(.+?){l}', text, re.DOTALL)[0].strip()
    m2 = re.sub(r'[^\S\r\n]{2,}', '\t', match)
    m3 = m2.replace('\n\t', '\n')
    df = pd.read_csv(StringIO(m3), sep="\t", names=columns, error_bad_lines=error_lines)
    return df
