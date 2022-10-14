def load_seq_scan(name, folder, haldane=False):
    '''
    load a seqscan result into a dataframe
    '''
    import pandas as pd
    pre = f'{folder}/{name}'
    info = pd.read_csv(f'{pre}.c100i100.tsv', sep='\t', names=['c100i100', 'pid'])
    info = info.merge(pd.read_csv(f'{pre}.c80i70.tsv', sep='\t', names=['c80i70', 'c100i100']), how="left")
    info = info.merge(pd.read_csv(f'{pre}.c80e3.tsv', sep='\t', names=['c80e3', 'c80i70']), how="left")
    if haldane:
        info = info.merge(pd.read_csv(f'{pre}.aravind.scan.arch', sep='\t', names=['c100i100', 'aravind'], usecols=[0,1], skiprows=[0]), how="left")
        info = info.merge(pd.read_csv(f'{pre}.pfam.hmmscan.arch', sep='\t', names=['c100i100', 'pfam'], usecols=[0,1], skiprows=[0]), how="left")
    else:    
        info = info.merge(pd.read_csv(f'{pre}.query.profiledb.rps.arch', sep='\t', names=['c100i100', 'aravind'], usecols=[0,1], skiprows=[0]), how="left")
        info = info.merge(pd.read_csv(f'{pre}.query.pfam.rps.arch', sep='\t', names=['c100i100', 'pfam'], usecols=[0,1], skiprows=[0]), how="left")
    return info


def cluster2aln(group_cluster,df,esl_index_file, grouper='c80e3', redundancy_cluster='c80i70',align_method='famsa', query=False, cpu=12):
    import os
    import tempfile
    from subprocess import Popen,PIPE
    from rotifer.devel.beta.sequence import sequence
    with tempfile.TemporaryDirectory() as tmpdirname:
        if query:
            df.query(query)[redundancy_cluster].drop_duplicates().dropna().to_csv(f'{tmpdirname}/accs', index=None, header=None)
        else:
            df[df[grouper]  == group_cluster][redundancy_cluster].drop_duplicates().dropna().to_csv(f'{tmpdirname}/accs', index=None, header=None)
        if not os.path.exists(esl_index_file+'.ssi'):
            Popen(f'esl-sfetch --index {esl_index_file}',stdout=PIPE, shell=True).communicate()
        Popen(f'esl-sfetch -f {esl_index_file} {tmpdirname}/accs > {tmpdirname}/accs.fa',stdout=PIPE, shell=True).communicate()
        b = sequence(f'{tmpdirname}/accs.fa').align(method=align_method, cpu=cpu)
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

def annotation(seqobj, coordinates, delimiter=True ):
    '''
    Method that recives an sequence object and a list with tuples containing [(start, end annotation)]
    To add annotaion in aligment
    Example:
        gf.annotation(seqobject, [(35, 217, ' domain 1 (probability 90%, WP_xxxxxx)'), (240, 380, 'domain 2')] )
    
    I sill have to finish it, but I  am almost done


    '''
    import pandas as pd
    s = seqobj.copy()
    t = pd.Series(list(s.df.iloc[0,1]))

    t.iloc[:] =' '
    for x in coordinates:
        if len(x) > 3:
            query_anchor = x[3]
            tq = pd.Series(list(s.df.query('id ==@query_anchor').sequence.values[0])) 
            start = tq.where(lambda x: x!='-').dropna().iloc[x[0]:x[1]].index.min() 
            end = tq.where(lambda x: x!='-').dropna().iloc[x[0]:x[1]].index.max() 
        else:
            start = x[0]
            end = x[1]

        annotation = x[2]
        size = end - start
        size_an = len(annotation)
        if size  -2 < size_an :            
                x =   size -size_an + -2  
                annotation = annotation[0:x]

        if delimiter:
             an = f'|{annotation}|'.center(size,'*')
        else:
            an = f'{annotation}'.center(size,'*')
        t.update(pd.Series(list(an), index=range(start,end)))
    s.df = pd.concat([pd.DataFrame([['annotation', ''.join(t.to_list()), 'annotation']], columns=['id', 'sequence', 'type']),s.df])
    return s

def hmmsearch_full2pandas (file, error_lines=True, keep_threshold=False ):
    """
    Function to load the result of hmmsearch or scan ang obtain the results for the full protein
    If error_lines = False, it skip the results that gave error loading to the dataframe.
    Those errors only happens if any field present an string with more than two spaces in sequence (mostly found in  descrption field )
    """
    import pandas as pd
    import re
    from io import StringIO
    
    with open(file, 'r') as f:
        text = f.read()
    columns = 'Evalue  score  bias    BD_Evalue  BD_score  BD_bias    exp  N  Sequence       Description'.split()
    #u = '------- ------ -----    ------- ------ -----   ---- --  --------       -----------'
    u =  '------- ------ -----    ------- ------ -----   ---- --  --------   -----------'
    l = 'Domain annotation for each sequence'
    match = re.findall(f'{u}(.+?){l}', text, re.DOTALL)[0].strip()
    m2 = re.sub(r'[^\S\r\n]{2,}', '\t', match)
    m3 = m2.replace('\n\t', '\n')
    df = pd.read_csv(StringIO(m3), sep="\t", names=columns,usecols = [0,1,2,3,4,5,6,7,8],  error_bad_lines=error_lines)
    if keep_threshold:
        return df
    
    df = df.query('~Evalue.str.contains("threshold")')
    df.Evalue = df.Evalue.astype(float)

    return df



def hhr_to_aln(seqobj, hhr, database=False):

    """
    Function to add hhr info into sequence object aligment
    """
    from rotifer.devel.alpha import gian_func as gf
    import rotifer.interval.utils as ru
    from collections import OrderedDict
    import re
    import pandas as pd


    si = re.findall(f'No ([0-9]+)', hhr, re.DOTALL)
    si = [int(x) for x in si]
    query = re.findall('Query\s+(.*?)\s',hhr, re.MULTILINE)[0]
    tdf = pd.DataFrame()
    for x in range(len(si)):
        try:
            match = re.findall(
                f'No {si[x]}(.+?)No ',
                hhr,
                re.DOTALL
            )[0].strip()
        except IndexError:
            match = re.findall(
                f'No {si[x]}(.+)\n',
                hhr,
                re.DOTALL
            )[0].strip()
        if database == "pfam":    
            try:
                idp = re.findall(
                    f'>.+?;(.+?);',
                    match,
                    re.DOTALL
                )[0].strip()
            except IndexError:
                idp = re.findall(
                    f'>(.+?)\n',
                    match,
                    re.DOTALL
                )[0].strip()
        else:        
            try:
                idp = re.findall(
                    f'>(.+?)\.',
                    match,
                    re.DOTALL
                )[0].strip()
            except IndexError:
                idp = re.findall(
                    f'>(.+?)\n',
                    match,
                    re.DOTALL
                )[0].strip()


        values = match.split('\n')[1].split()
        values =pd.Series(values).str.split(
            '=',
            expand=True
        ).rename({
            0: 'col',
            1:'val'},
                 axis=1)
        values['model'] = idp
        values = values.pivot(
            index='model',
            columns ='col',
            values = 'val'
        )
        values['query'] = query
        values['start'] = re.findall(
            f'{query}\s+([0-9]+) ',
            match,
            re.DOTALL
        )[0]
        values['end'] = re.findall(
            f'{query}.+?([0-9]+?) \(',
            match,
            re.DOTALL)[-1]

        tdf = pd.concat([tdf, values])

    hhr_nolr = ru.filter_nonoverlapping_regions(
        tdf.astype({'start':int, 'end':int}).reset_index(),
        end='end',
        start='start',
        reference='query',
        criteria=OrderedDict([('Probab', False), ('region_length', True)]
                             ))
    # ISSUE !!!!!! When the rep seq do not have the first domains
    # it breaks the code because there is no way to map the domain.
    # This is an imporvisation that will overpass it, but with error in annotation!!
    hhr_nolr = hhr_nolr[hhr_nolr.end > hhr_nolr.start]
    return hhr_nolr



def add_arch_to_df(df):
    '''
    Add architecture and clusters info to a neighborhood df
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    df = df.copy() 
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        sequence(
                df.query(
                    'type =="CDS"'
                    ).pid.dropna().drop_duplicates().to_list()
                ).to_file(f'{tmpdirname}/seqfile') 
        

        os.chdir(tmpdirname)
        os.makedirs('tmpd')
        Popen(f'/home/nicastrogg/bin/scanseqs.nih.sh tmp seqfile' , stdout=PIPE,shell=True).communicate()
        info = pd.read_csv(f'tmp.c100i100.tsv', sep='\t', names=['c100i100', 'pid'])
        info = info.merge(pd.read_csv(f'tmp.c80i70.tsv', sep='\t', names=['c80i70', 'c100i100']), how="left")
        info = info.merge(pd.read_csv(f'tmp.c80e3.tsv', sep='\t', names=['c80e3', 'c80i70']), how="left")
        info = info.merge(pd.read_csv(f'tmp.query.profiledb.rps.arch', sep='\t', names=['c100i100', 'profiledb'], usecols=[0,1], skiprows=[0]), how="left")
        info = info.merge(pd.read_csv(f'tmp.query.pfam.rps.arch', sep='\t', names=['c100i100', 'pfam'], usecols=[0,1], skiprows=[0]), how="left")
        df = df.merge(info, how='left')
        os.chdir(cwd)
        
    return df


def annotate_seqobj(seqobj,df,cnt='profiledb'):
    '''
    Add annotation to seqobject tax,clusters,compact_neighborhood 
    '''
    df = df.copy() 
    seqobj = seqobj.copy()
    accs = seqobj.df.id.to_list()
    cn = df[df.block_id.isin(df.query('pid in @accs').block_id)].compact_neighborhood(cnt)
    seqobj.df = seqobj.df.merge(df.query('pid in @accs')[['pid','block_id','profiledb', 'pfam', 'c80e3', 'c80i70']].rename({'pid':'id'},axis=1).set_index('block_id').join(cn), on='id', how='left')

    return seqobj



def psiblast(acc,db='nr.50', cpu=96):
    '''
    Psiblast it can accept sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cols =['hit','query','hitstart','hitend','evalue','querycoverage','querystart','queryend','iteration', 'bitscore' ,'length']
    if db.startswith('nr'):
        db = f'{os.environ["FADB"]}/nr/{db}.mmseqs.fa'
    else:
        db = f'{os.environ["FADB"]}/all.fa'

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        if isinstance (acc, sequence):
            acc.to_file(f'{tmpdirname}/seqfile') 
            Popen(f'splishpsi -a {cpu} -in_msa {tmpdirname}/seqfile -d {db} > {tmpdirname}/out' , stdout=PIPE,shell=True).communicate()
            Popen(f'blast2table {tmpdirname}/out > {tmpdirname}/out.tsv' , stdout=PIPE,shell=True).communicate()
            t = pd.read_csv(f'{tmpdirname}/out.tsv', sep='\t', names=cols)
            with open(f'{tmpdirname}/out') as f:
                    blast_r = f.read()
        
    return (t, blast_r) 



def search2aln(df, coverage=50, evalue=1e-3, arch=None):
    import os
    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import pandas as pd
    if coverage > 1:
        coverage = coverage/100
    '''
    From a search tabular file, filter by cov, and evalue to make a aln object
    '''
    df = df.query(f'evalue <= {evalue} and querycoverage > {coverage}')
    seqobj = sequence(df.hit.drop_duplicates().to_list())
    if arch == 'profiledb':
        arch = ' '
    if arch:
        with tempfile.TemporaryDirectory() as tmpdirname:
            seqobj.to_file(f'{tmpdirname}/fasta')
            Popen(f'cat {tmpdirname}/fasta|splishrps {arch}|rps2arch -pcut 0.001 -overlap 0.03 > {tmpdirname}/arch', stdout=PIPE,shell=True).communicate()
            archdf = pd.read_csv(f'{tmpdirname}/arch',sep="\t", names=['id', 'arch', 'position'])

            

    seqobj.df = seqobj.df.merge(df[['hit','evalue','hitstart','hitend']].rename({'hit':'id'},axis=1), how='left')
    seqobj.df.sequence = seqobj.df.apply(lambda x: x.sequence[x.hitstart:x.hitend], axis=1)
    seqobj = seqobj.align()
    seqobj.df = seqobj.df.merge(df[['hit','evalue','hitstart','hitend']].rename({'hit':'id'},axis=1), how='left')
    seqobj = seqobj.sort(['evalue'])
    if arch:
        seqobj.df = seqobj.df.merge(archdf, how='left')

    seqobj.df.id = seqobj.df.apply(lambda x: f'{x.id}/{x.hitstart}-{x.hitend}', axis=1)

    return seqobj


def add_arch_to_seqobj(seqobj,db='profiledb', cpu=96):
    '''
    Psiblast it can accept sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cols =['pid','arch','evalue']
    if db == 'profiledb':
        db = ' '

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        acc = sequence(seqobj.df.id.to_list())
        acc.to_file(f'{tmpdirname}/seqfile') 
        Popen(f'cat {tmpdirname}/seqfile| splishrps -a {cpu} {db} > {tmpdirname}/out' , stdout=PIPE,shell=True).communicate()
        Popen(f'rps2arch{tmpdirname}/out > {tmpdirname}/out.tsv' , stdout=PIPE,shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/out.tsv', sep='\t', names=cols)
        seqobj.df = seqobj.df.merge(t, how='left) 
    return seqobj
