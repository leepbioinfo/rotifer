def load_seq_scan(name, folder, haldane=False):
    '''
    load a seqscan result into a dataframe
    '''
    import pandas as pd
    pre = f'{folder}/{name}'
    info = pd.read_csv(
            f'{pre}.c100i100.tsv', sep='\t', names=['c100i100', 'pid']
            )
    info = info.merge(
            pd.read_csv(
                f'{pre}.c80i70.tsv',
                sep='\t',
                names=['c80i70', 'c100i100']
                ), how="left"
            )
    info = info.merge(
            pd.read_csv(
                f'{pre}.c80e3.tsv',
                sep='\t',
                names=['c80e3', 'c80i70']
                ),
            how="left")
    if haldane:
        info = info.merge(
                pd.read_csv(
                    f'{pre}.aravind.scan.arch',
                    sep='\t',
                    names=['c100i100', 'aravind'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
        info = info.merge(
                pd.read_csv(
                    f'{pre}.pfam.hmmscan.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
    else:
        info = info.merge(
                pd.read_csv(
                    f'{pre}.query.profiledb.rps.arch',
                    sep='\t',
                    names=['c100i100', 'aravind'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
        info = info.merge(
                pd.read_csv(
                    f'{pre}.query.pfam.rps.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
    return info


def cluster2aln(
        group_cluster,
        df,
        esl_index_file,
        grouper='c80e3',
        redundancy_cluster='c80i70',
        align_method='famsa',
        query=False,
        cpu=12
        ):
    import os
    import tempfile
    from subprocess import Popen, PIPE
    from rotifer.devel.beta.sequence import sequence
    with tempfile.TemporaryDirectory() as tmpdirname:
        if query:
            df.query(
                    query
                    )[redundancy_cluster].drop_duplicates().dropna().to_csv(
                            f'{tmpdirname}/accs',
                            index=None,
                            header=None
                            )
        else:
            df[df[grouper] == group_cluster][
                    redundancy_cluster
                    ].drop_duplicates(
                            ).dropna(
                                    ).to_csv(
                                            f'{tmpdirname}/accs',
                                            index=None,
                                            header=None
                                            )
        if not os.path.exists(esl_index_file+'.ssi'):
            Popen(
                    f'esl-sfetch --index {esl_index_file}',
                    stdout=PIPE,
                    shell=True
                    ).communicate()
        Popen(
                f'esl-sfetch -f {esl_index_file} {tmpdirname}/accs > \
                        {tmpdirname}/accs.fa',
                stdout=PIPE,
                shell=True
                ).communicate()
        b = sequence(
                f'{tmpdirname}/accs.fa'
                ).align(
                        method=align_method,
                        cpu=cpu)
        return b


def chunks(list_to_chunk, chunk_size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(list_to_chunk), chunk_size):
        yield list_to_chunk[i:i + chunk_size]


def cluster_Co_occurrence(
        df,
        count='c80e3',
        freq_cutoff=0.3,
        only_query=True,
        annotation=False):
    ''' Function to count the co-occurence of clusters within
    the NeighborhoodDF.
    The count parameter shoul use to define the cluster
    that would be analysed and
    the freq_cutoff is to use to define a minimum cut_off
    to display
    '''
    x = df[['block_id', count]].merge(
            df[['block_id', count]],
            how='outer',
            on='block_id'
            )
    x.columns = [
            'block_id',
            'query_cluster',
            'neighbor_cluster'
            ]
    x = x[x.query_cluster != x.neighbor_cluster]
    if only_query:
        xx = x[
                x['query_cluster'].isin(
                    df.query('query ==1').pid.unique()
                    )
                ]
    else:
        xx = x

    xxx = xx.groupby(
            ['query_cluster', 'neighbor_cluster']
            ).block_id.nunique().reset_index()
    xxxx = xx.groupby(
            'query_cluster'
            ).block_id.nunique().rename(
                    'query_blocks'
                    ).reset_index().merge(
                            xxx,
                            how='left'
                            ).sort_values(
                                    ['query_blocks', 'block_id'],
                                    ascending=False
                                    )
    xxxx['query_freq'] = xxxx.block_id/xxxx.query_blocks
    if not annotation:
        return xxxx.query('query_freq >= @freq_cutoff')

    andf = df.groupby(count).agg(
            pfam=('pfam', count_series),
            aravind=('aravind', count_series)
            ).reset_index()
    xxxx = xxxx.merge(
            andf.rename({
                count: 'query_cluster',
                'pfam': 'query_pfam',
                'aravind': 'query_aravind'
                }, axis=1),
            how='left')
    xxxx = xxxx.merge(
            andf.rename({
                count: 'neighbor_cluster',
                'pfam': 'neighbor_pfam',
                'aravind': 'neighbor_aravind'
                }, axis=1),
            how='left')
    return xxxx.query('query_freq >= @freq_cutoff')


def count_series(
        series,
        normalize=False,
        cut_off=False,
        count='domain'
        ):
    '''
    Function to flatten a pd.Series.value_counts
    and count the number of domains recognized in
    the whole proteins set of proteins
    Mind that by default, it count domain not architecture.
    If one wants to count architecture, it should change the count
    option to 'architecture'.
    The normalize options is to write the results as frequency.
    The cut_off options is to print results only above a given threashold.
    If the normalize function is True, you should give the cutoff value as pct.
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
    flattened_list = []
    if count == 'architecture':
        s = series.dropna().value_counts(normalize=normalize)
    else:
        s = series.str.split('+').explode().dropna().value_counts(
                normalize=normalize
                )

    if cut_off:
        cut_off = cut_off/100
        s = s.where(lambda x: x >= cut_off).dropna()
    for y, z in s.items():
        if normalize:
            flattened_list.append(f'{y}({100 * z:.2f}%)')
        else:
            flattened_list.append(f'{y}({z})')
    return ', '.join(flattened_list)


def fetch_seq(seqs):
    import time
    from Bio import Entrez
    from rotifer.db.ncbi import NcbiConfig
    from rotifer.devel.beta.sequence import sequence
    from rotifer.devel.alpha.gian_func import chunks
    if isinstance(seqs, list):
        seqs = chunks(seqs, 200)
    else:
        return print('add seuquences as list object')

    seq_string = ''
    for x in seqs:
        f = Entrez.efetch(
                db='protein',
                rettype='fasta',
                retmode='text',
                id=','.join(x),
                api_key=NcbiConfig['api_key']
                ).read()
        seq_string = seq_string + f
        time.sleep(1)
    return sequence.from_string(seq_string)


def annotation(seqobj, coordinates, delimiter=True):
    '''
    Method that receives an sequence object and a list with tuples
    containing [(start, end annotation)]
    To add annotaion in aligment
    Example:
        gf.annotation(
        seqobject,
        [(35, 217,
        ' domain 1 (probability 90%, WP_xxxxxx)'),
        (240, 380, 'domain 2')
        ] )
    I sill have to finish it, but I  am almost done
    '''
    import pandas as pd
    s = seqobj.copy()
    t = pd.Series(list(s.df.iloc[0, 1]))

    t.iloc[:] = ' '
    for x in coordinates:
        if len(x) > 3:
            query_anchor = x[3]
            tq = pd.Series(
                    list(
                        s.df.query(
                            'id ==@query_anchor'
                            ).sequence.values[0]
                        )
                    )
            start = tq.where(
                    lambda x: x != '-'
                    ).dropna().iloc[x[0]:x[1]].index.min()
            end = tq.where(
                    lambda x: x != '-'
                    ).dropna().iloc[x[0]:x[1]].index.max()
        else:
            start = x[0]
            end = x[1]

        annotation = x[2]
        size = end - start
        size_an = len(annotation)
        if (size - 2) < size_an:
            x = size - size_an + - 2
            annotation = annotation[0:x]

        if delimiter:
            an = f'|{annotation}|'.center(size, '*')
        else:
            an = f'{annotation}'.center(size, '*')
        t.update(pd.Series(list(an), index=range(start, end)))
    s.df = pd.concat(
            [pd.DataFrame(
                [['annotation', ''.join(t.to_list()), 'annotation']],
                columns=['id', 'sequence', 'type']
                ),
             s.df
             ])
    return s


def hmmsearch_full2pandas(
        file,
        error_lines=True,
        keep_threshold=False
       ):
    """
    Function to load the result of hmmsearch or scan ang obtain
    the results for the full protein.
    If error_lines = False, it skip the results that gave error
    loading to the dataframe.
    Those errors only happens if any field present an string with
    more than two spaces in sequence (mostly found in  descrption field )
    """
    import pandas as pd
    import re
    from io import StringIO
    with open(file, 'r') as f:
        text = f.read()
    columns = ['Evalue',
               'score',
               'bias',
               'BD_Evalue',
               'BD_score',
               'BD_bias',
               'exp',
               'N',
               'Sequence',
               'Description'
               ]
    # u = '------- ------ -----    ------- ------ -----   ---- --  --------       -----------'
    u = '------- ------ -----    ------- ------ -----   ---- --  --------   -----------'
    string_tore = 'Domain annotation for each sequence'
    match = re.findall(f'{u}(.+?){string_tore}', text, re.DOTALL)[0].strip()
    m2 = re.sub(r'[^\S\r\n]{2,}', '\t', match)
    m3 = m2.replace('\n\t', '\n')
    df = pd.read_csv(
            StringIO(m3),
            sep="\t",
            names=columns,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            error_bad_lines=error_lines)
    if keep_threshold:
        return df

    df = df.query('~Evalue.str.contains("threshold")')
    df.Evalue = df.Evalue.astype(float)

    return df


def hhr_to_aln(seqobj, hhr, database=False):

    """
    Function to add hhr info into sequence object aligment
    """
    import rotifer.interval.utils as ru
    from collections import OrderedDict
    import re
    import pandas as pd

    si = re.findall('No ([0-9]+)', hhr, re.DOTALL)
    si = [int(x) for x in si]
    query = re.findall(r'Query\s+(.*?)\s', hhr, re.MULTILINE)[0]
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
                    r'>.+?;(.+?);',
                    match,
                    re.DOTALL
                )[0].strip()
            except IndexError:
                idp = re.findall(
                    r'>(.+?)\n',
                    match,
                    re.DOTALL
                )[0].strip()
        else:
            try:
                idp = re.findall(
                    r'>(.+?)\.',
                    match,
                    re.DOTALL
                )[0].strip()
            except IndexError:
                idp = re.findall(
                    r'>(.+?)\n',
                    match,
                    re.DOTALL
                )[0].strip()

        values = match.split('\n')[1].split()
        values = pd.Series(values).str.split(
            '=',
            expand=True
        ).rename({
            0: 'col',
            1: 'val'},
                 axis=1)
        values['model'] = idp
        values = values.pivot(
            index='model',
            columns='col',
            values='val'
        )
        values['query'] = query
        values['start'] = re.findall(
            rf'{query}\s+([0-9]+) ',
            match,
            re.DOTALL
        )[0]
        values['end'] = re.findall(
            rf'{query}.+?([0-9]+?) \(',
            match,
            re.DOTALL)[-1]

        tdf = pd.concat([tdf, values])

    hhr_nolr = ru.filter_nonoverlapping_regions(
        tdf.astype({'start': int, 'end': int}).reset_index(),
        end='end',
        start='start',
        reference='query',
        criteria=OrderedDict([('Probab', False), ('region_length', True)]
                             ))
    # ISSUE !!!!!! When the rep seq do not have the first domains
    # it breaks the code because there is no way to map the domain.
    # This is an imporvisation that will overpass it,
    # but with error in annotation!!
    hhr_nolr = hhr_nolr[hhr_nolr.end > hhr_nolr.start]
    return hhr_nolr


def add_arch_to_df(seqobj, full=False, inplace=True):
    '''
    Add architecture and clusters info to a neighborhood df
    '''

    import tempfile
    from subprocess import Popen, PIPE
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    seqobjc = seqobj.copy()
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdirname:
        if hasattr(seqobjc, 'ndf'):
            df = seqobjc.ndf
        # temporary save fasta sequence file
            sequence(
                    df.query(
                        'type =="CDS"'
                        ).pid.dropna().drop_duplicates().to_list()
                    ).to_file(f'{tmpdirname}/seqfile')
        else:
            print('Seq object do not have neighborhood df, running annotation on aligment sequences')
            df = seqobjc.df
            sequence(
                    df.id.dropna().drop_duplicates().to_list()
                    ).to_file(f'{tmpdirname}/seqfile')

        os.chdir(tmpdirname)
        os.makedirs('tmpd')
        if full:
            script = 'scanseqs.nih.sh'
        else:
            script = 'scanseqs.short.nih.sh' 

        Popen(
                f'/home/nicastrogg/bin/{script} tmp seqfile',
                stdout=PIPE,
                shell=True
                ).communicate()
        info = pd.read_csv(
                'tmp.c100i100.tsv',
                sep='\t',
                names=['c100i100', 'pid']
                )
        info = info.merge(
                pd.read_csv(
                    'tmp.c80i70.tsv',
                    sep='\t',
                    names=['c80i70', 'c100i100']
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.c80e3.tsv',
                    sep='\t',
                    names=['c80e3', 'c80i70']),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.query.profiledb.rps.arch',
                    sep='\t',
                    names=['c100i100', 'profiledb'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.query.pfam.rps.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.query.both.rps.arch',
                    sep='\t',
                    names=['c100i100', 'arch'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        if inplace:
            if hasattr(seqobjc, 'ndf'):
                seqobj.ndf = df.merge(info, how='left')
                seqobj.ndf = seqobj.ndf.drop_duplicates().reset_index(drop=True)
            else:
                info = info.rename({'pid':'id'}, axis=1)
                seqobj.df = df.merge(info, how='left')

            seqobj.profiledb = pd.read_csv(
                    'tmp.query.profiledb.rps.dom.tsv',
                    sep='\t', names=['pid', 'model', 'star', 'end'])
            with open('tmp.profiledb.rps.query.out') as f:
                seqobj.profiledbout = ''.join(f.readlines())
            seqobj.pfam = pd.read_csv(
                    'tmp.query.pfam.rps.dom.tsv',
                    sep='\t', names=['pid', 'model', 'start', 'end'])
            with open('tmp.pfam.rps.query.out') as f:
                seqobj.pfamout = ''.join(f.readlines())
            seqobj.arch = pd.read_csv(
                    'tmp.query.both.rps.dom.tsv',
                    sep='\t', names=['pid', 'model', 'start', 'end'])
            os.chdir(cwd)
            return 'architecuture add to seqobj'

        ### The annotation of arch on sequence objct is not working if not inplace ########
        with open('tmp.pfam.rps.query.out') as f:
            seqobjc.pfamout = ''.join(f.readlines())

        seqobjc.pfam = pd.read_csv(
                'tmp.query.pfam.rps.dom.tsv',
                sep='\t')

        seqobjc.arch = pd.read_csv(
                'tmp.query.both.rps.dom.tsv',
                sep='\t')

        seqobjc.profiledb = pd.read_csv(
                'tmp.query.profiledb.rps.dom.tsv',
                sep='\t')
        with open('tmp.profiledb.rps.query.out') as f:
            seqobjc.profiledbout = ''.join(f.readlines())
        seqobjc.ndf = df.merge(info, how='left')
        seqobjc.ndf = seqobj.ndf.drop_duplicates().reset_index(drop=True)
        os.chdir(cwd)
        
    return seqobjc


def annotate_seqobj(seqobj,
                    ndf = False,
                    cnt='arch',
                    selected_collumns=False
                    ):
    '''
    Add annotation to seqobject tax,clusters,compact_neighborhood
    '''

    import pandas as pd

    if isinstance(ndf, pd.DataFrame):
        df = ndf.copy()
    else:
        df = seqobj.ndf.copy()
    seqobj = seqobj.copy()
    accs = seqobj.df.id.to_list()
    cn = df[
            df.block_id.isin(df.query('pid in @accs').block_id)
            ].compact_neighborhood(cnt)
    if selected_collumns:
        selected_collumns = selected_collumns
    else:
        selected_collumns = ['pid',
         'block_id',
         'profiledb',
         'pfam',
         'arch',                    
         'c80e3',
         'c80i70']

    seqobj.df = seqobj.df.merge(df.query('pid in @accs')[selected_collumns ].rename(
            {'pid': 'id'}, axis=1
            ).set_index('block_id').join(cn), on='id', how='left')

    return seqobj


def psiblast(acc,
             db='nr.50',
             cpu=96,
             aln=True,
             max_out = 5000):
    '''
    Psiblast it can accept sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cols =[
            'hit',
            'query',
            'hitstart',
            'hitend',
            'evalue',
            'querycoverage',
            'querystart',
            'queryend',
            'iteration',
            'bitscore',
            'length'
            ]
    print(db)
    if db.startswith('nr'):
        db = f'{os.environ["FADB"]}/nr/{db}.mmseqs.fa'
    elif db.startswith('all'):
        db = f'{os.environ["FADB"]}/allfa/{db}.mmseqs.fa'
    elif db == "prok" :
        db = f'{os.environ["FADB"]}/prok.fa '
    elif db == "euk":
        db = f'{os.environ["FADB"]}/euk.fa'
    elif db == "af":
        db = f'{os.environ["FADB"]}/uniref50.fa'
    else:
        db = f'{os.environ["FADB"]}/all.fa'

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        if isinstance (acc, sequence):
            acc.to_file(f'{tmpdirname}/seqfile') 
            if aln:
                Popen(f'splishpsi -a {cpu} -in_msa {tmpdirname}/seqfile -d {db} -num_alignments {max_out} -num_descriptions {max_out} > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            else:
                Popen(f'splishpsi -a {cpu} -i {tmpdirname}/seqfile -d {db} -num_alignments {max_out} -num_descriptions {max_out} > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            Popen(f'blast2table {tmpdirname}/out > {tmpdirname}/out.tsv',
                  stdout=PIPE,
                  shell=True).communicate()
            t = pd.read_csv(f'{tmpdirname}/out.tsv', sep="\t", names=cols)
            with open(f'{tmpdirname}/out') as f:
                    blast_r = f.read()
        
    return (t, blast_r) 



def search2aln(df, coverage=50, evalue=1e-3, id_with_coord = False, arch=None):
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
    seqobj.df.id = seqobj.df.apply(lambda x: f'{x.id}/{x.hitstart}-{x.hitend}', axis=1)
    evalues = seqobj.df.copy()[['id', 'evalue']]
    seqobj = seqobj.align()
    seqobj.df = seqobj.df.merge(evalues,on=['id'], how='left')
    seqobj = seqobj.sort(['evalue'])
    if arch:
        # We should fix the arch with coordinates similar to the evalue strategy
        print('arch option may contains bug!i')
        seqobj.df = seqobj.df.merge(archdf, how='left')

    #seqobj.df.id = seqobj.df.apply(lambda x: f'{x.id}/{x.hitstart}-{x.hitend}', axis=1)

    if id_with_coord:
        seqobj.df = seqobj.df.drop_duplicates(['id']).reset_index(drop=True)
        return seqobj
    seqobj.df[['hit_start', 'hit_end']] = seqobj.df.id.str.split('/', expand=True)[1].str.split('-', expand=True)
    seqobj.df.id = seqobj.df.id.str.split('/', expand=True)[0]
    seqobj.df = seqobj.df.drop_duplicates(['id', 'hit_start', 'hit_end']).reset_index(drop=True)
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
    cols =['id','arch','evalue']
    if db == 'profiledb':
        db = ' '

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        acc = sequence(seqobj.df.id.to_list())
        acc.to_file(f'{tmpdirname}/seqfile') 
        Popen(f'cat {tmpdirname}/seqfile| splishrps -a {cpu} {db} > {tmpdirname}/out' , stdout=PIPE,shell=True).communicate()
        Popen(f'rps2arch {tmpdirname}/out > {tmpdirname}/out.tsv' , stdout=PIPE,shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/out.tsv', sep='\t', names=cols)
        seqobj.df = seqobj.df.merge(t, how='left', on='id') 
    return seqobj


def full_annotate(seqobj,
                  progress=True,
                  batch_size=8,
                  mirror="/am/ftp-genomes",
                  threads=None,
                  after=5,
                  before=5,
                  eukaryotes=False,
                  full=True):
    from rotifer.db import ncbi 
    from rotifer.devel.alpha import gian_func as gf
    gnc = ncbi.GeneNeighborhoodCursor(
            progress=progress,
            batch_size=batch_size,
            mirror=mirror,
            threads=threads,
            after=after,
            before=before,
            eukaryotes=eukaryotes)
    seqobj.ndf = gnc.fetchall(seqobj.df.id.to_list())
    gf.add_arch_to_df(seqobj, full=full)
    return seqobj

def padding_df(df, how='right'):
    cdf = df.copy()
    c = df.columns
    pad_col_name=[]
    for x in c:
        cdf[x] = cdf[x].fillna('').astype(str)
        w = cdf[x].str.len().max()
        cdf[x] = cdf[x].str.pad(width =w, side=how)
        if how == 'left':
            pad_col_name.append(x.rjust(w))
        else:
            pad_col_name.append(x.ljust(w))

    cdf.columns = pad_col_name    
    return cdf


def alnxaln(seqobj, clustercol = 'c50i0', minseq=10):
    import os
    import tempfile
    import pandas as pd
    from rotifer.io import hhsuite
    
    path = os.getcwd()
    result_table = pd.DataFrame()
    if clustercol in seqobj.df.columns:
        pass
    else:
       print(f'add the cluster{clustercol} columns to your seqobj')
       return

    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        seqobj = seqobj.copy()
        alndict = dict()
        clist = list(seqobj.df[clustercol].value_counts().where(lambda x: x>= minseq).dropna().index)
        for x in clist:
            alndict[x] = seqobj.filter(f'{clustercol} ==\"{x}\"').align(method='linsi')
            alndict[x].to_file(f'{x}.aln') 
        for x in clist:
            os.system(f'for x in *.aln; do hhalign -i $x -t {x}.aln;done')
            os.system('cat *.hhr >> allhhr.txt')
            hhrs = hhsuite.read_hhr('./')
            result_table = pd.concat([result_table, hhrs])
        with open('./allhhr.txt', 'r') as f : allhhr = f.read()        
    os.chdir(path)
    return (result_table, allhhr, alndict) 



def split_by_model(seqobj, model):
    '''
    Using a hmm model to split your sequence object to match only the model region.
    If more than one match in one protein, it will split the match in two sequence.
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    import rotifer
    from rotifer.core import loadpath
    
    hmmer2table_from_rotifer = f'{loadpath.path().globalpath}rotifer/bin/hmmer2table'

    nseqobj =seqobj.copy()
    with tempfile.TemporaryDirectory() as tmpdirname:
        nseqobj.df.sequence = nseqobj.df.sequence.str.replace('-','')
        nseqobj.to_file(f'{tmpdirname}/seqfile') 

        if isinstance(model, rotifer.devel.beta.sequence.sequence):
            model.to_file(f'{tmpdirname}/model')
            Popen(f'hmmbuild {tmpdirname}/model.hmm {tmpdirname}/model',
                  stdout=PIPE,
                  shell=True).communicate()
            model = f'{tmpdirname}/model.hmm'

        Popen(f'hmmsearch {model} {tmpdirname}/seqfile > {tmpdirname}/hhsearch_r',
              stdout=PIPE,
              shell=True).communicate()
        Popen(f'{hmmer2table_from_rotifer} {tmpdirname}/hhsearch_r |domain2architecture |architecture2table > {tmpdirname}/hhsearch_table',
              stdout=PIPE,
              shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/hhsearch_table', sep="\t")

    nseqobj.df = nseqobj.df.merge(t.rename({'ID':'id'},axis=1), how='left')
    nseqobj.df.end = nseqobj.df.end.fillna(nseqobj.df.length)
    nseqobj.df.start = nseqobj.df.start.fillna(1)
    nseqobj.df.sequence = nseqobj.df.apply(lambda x: x.sequence[int(x.start):int(x.end)], axis=1)
    nseqobj = nseqobj.filter('domain =="model"').align()
    return (nseqobj)

def remove_redundancy(seqobj, identity =80, coverage = 70):
    from rotifer.devel.beta.sequence import sequence as sequence
    import pandas as pd
    s = seqobj.copy()
    if len(s.df.id.value_counts().value_counts()) >1:
        print('multiples sequences with same name, using index to remove redundancy')

    s.add_cluster(coverage=coverage, identity=identity, inplace=True)
    s = s.filter(keep=s.df[f'c{coverage}i{identity}'].drop_duplicates().to_list())
    return(s)

def add_cordinates_to_aln(seqobj):
    from rotifer.devel.beta.sequence import sequence as sequence
    c = seqobj.copy()
    cx = sequence(c.df.id.tolist())
    cx. df = cx.df.rename({'sequence':'full_sequence', 'length':'full_length'}, axis=1)
    c.df = c.df.merge(cx.df[['id','full_sequence','full_length']], how='left')
    c.df['start'] = c.df.apply(lambda x : 0 + x.full_sequence.find(x.sequence.replace('-','')), axis=1)
    c.df['end'] = c.df.start + c.df['length']
    c.df['C_term'] = c.df['full_length'] - c.df['end']
    c.df =  c.df.drop(['full_sequence'], axis=1)
    return c

def trim_unk_neigh(df, ann='profiledb'):
    '''
    Using a given column, it will trim the neighborhoood to do not have un annotated protein at the boarders.
    ann paramether is the collumn to check the unknow limits
    '''

    def xxx (idf):
        imin = idf.query(f'{ann} != "?"').index[0]
        imax = idf.query(f'{ann} != "?"').index[-1]
        return idf.loc[imin:imax]
    int_df = df.groupby(['block_id']).apply(xxx).reset_index(drop=True)


    return int_df

def extend_aln(seqobj, n_terminal=50, c_terminal=50):
    from rotifer.devel.beta import sequence as rdbs
    from rotifer.devel.alpha import gian_func as rdagf
    import numpy as np
    """
    Extend the N and/or C termianl parts of a given sequence alignment
    """
    tm = rdagf.add_cordinates_to_aln(seqobj)
    tm2 = rdbs.sequence(seqobj.df.id.tolist())
    tm2.df = tm2.df.merge(tm.df[['id', 'full_length', 'start', 'end', 'C_term']], how='left')
    tm2.df.end = np.where(tm2.df.end + c_terminal <= tm2.df.full_length, tm2.df.end + c_terminal, tm2.df.full_length)
    tm2.df.start = np.where(tm2.df.start - n_terminal > 0, tm2.df.start - n_terminal, 1)
    tm2.df.sequence = tm2.df.apply(lambda x : x.sequence[x.start: x.end], axis=1)
    return tm2.align()

def get_correspondent_position(seqobj_source,seqobj_target, pid, position):
    """ Get the correspondent amminoacid position from one aligment/single sequence (rotifer seqobj) in a different sequence aligment
    seqobj_source: Sequence object where the user knows the desired position
    seqobt_target: Sequence object where the user need to find the correspondent position
    pid: Sequence used as anchor to find the correspondent position, it should be present in both source and target sequence object
    position: The postion in the pid wherethe user wants to find the correspondet position in the target sequence objct
    """
    from Bio import pairwise2
    import pandas as pd

    source =  seqobj_source.df.query('id == @pid').sequence.iloc[0].replace('-', '')
    target =  seqobj_target.df.query('id == @pid').sequence.iloc[0].replace('-','')
    aln = pairwise2.align.globalxx(source, target)[0]
    source = aln[0]
    target = aln[1]
    d = {'seq': list(source), 'trimmed' :list(target)}
    aln_df = pd.DataFrame(data=d).reset_index().rename({'index': 'pos'}, axis=1).query('trimmed != "-"').reset_index(drop=True)
    return int(aln_df.query('pos ==@position -1').index[0]) +1


def read_predicted_topologies(file):
    from rotifer.devel.beta import sequence as rdbs
    import pandas as pd
    import numpy as np
    def g (s):
        tomap = pd.Series(list(s.sequence))
        o = tomap.where(lambda x: x!="-").dropna().rename('aln').to_frame()
        o['top'] = list(s.topology)
        od = o.top.to_dict()
        tomap.update(od)
        return "".join(tomap.tolist())

    a = np.array(['header', 'sequence', 'topology'])
    df = pd.read_csv(file, names=['data'])
    df['data_type'] = np.resize(a, len(df))
    h = df.query('data_type =="header"').data.str.strip('>').rename("id").reset_index(drop=True)
    s = df.query('data_type =="sequence"').data.rename("sequence").reset_index(drop=True)
    t = df.query('data_type =="topology"').data.rename("topology").reset_index(drop=True)
    df2 = pd.concat([h,s,t], axis=1)
    aln = rdbs.sequence(df2[['id', 'sequence']]).align()
    aln.df['topology'] = df2.topology
    aln.df['aln_top'] = aln.df.apply(g, axis=1)
    aln.con = aln.add_consensus().df.iloc[0:4,0:2].rename({'sequence': 'aln_top'}, axis=1)
    aln.topology = pd.concat([aln.con,aln.df[['id', 'aln_top']]]).reset_index(drop=True)
    return aln

def uniref50_add_info(seqobj):
    seqobj = seqobj.copy()
    from rotifer.devel.beta.sequence import sequence
    from rotifer.db.local import ete3
    from rotifer.devel.alpha import gian_func as gf
    dd = sequence(seqobj.df.id.tolist())
    dd.add_cluster(coverage=0, identity=0, inplace=True)
    seqobj.df = seqobj.df.merge( dd.df[['id', 'length', 'c0i0']].rename({'length': 'plen'},axis=1))
    xx = ete3.TaxonomyCursor()
    seqobj.df[['description', 'cluster_size', 'Tax', 'TaxID', 'RepID']] = seqobj.df.description.str.split("n=|Tax=|TaxID=|RepID=", expand=True)
    seqobj = gf.add_cordinates_to_aln(seqobj)
    tid = xx.fetchall(seqobj.df.TaxID.drop_duplicates().tolist())
    tid.taxid = tid.taxid.astype(int)
    tid.rename({'taxid':'TaxID'}, axis=1, inplace=True)
    seqobj.df.TaxID = seqobj.df.TaxID.astype(int)
    seqobj.df = seqobj.df.merge(tid, how='left')
    return seqobj


def uniref50_to_ncbi(uniref50_strings):
     import pandas as pd
     import sqlite3
    # Connect to the SQLite database
     conn = sqlite3.connect('/netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db')
     cursor = conn.cursor()
     pids = [uniprotid.split('_')[1] for uniprotid in uniref50_strings]

     placeholders = ', '.join(['?' for _ in pids])
     sql_query = f"""
         SELECT *
         FROM EMBL
         WHERE uniprotid IN ({placeholders});
     """
     uniref50DF = pd.read_sql_query(sql_query, conn, params=pids)
     uniref50DF['Uniref50'] = 'UniRef50_' + uniref50DF.uniprotid
     # Close the cursor and connection
     cursor.close()
     conn.close()

     # Return the result as a pandas Series
     return uniref50DF

def uniref50_to_clusters(uniref50_strings):
     import pandas as pd
     import sqlite3
    # Connect to the SQLite database
     conn = sqlite3.connect('/netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db')
     cursor = conn.cursor()

     placeholders = ', '.join(['?' for _ in uniref50_strings])
     sql_query = f"""
         SELECT *
         FROM uniref50
         WHERE uniref50 IN ({placeholders});
     """
     uniref50DF = pd.read_sql_query(sql_query, conn, params=uniref50_strings)
     uniprotid = uniref50DF.uniprotid.unique().tolist() 

     placeholders2 = ', '.join(['?' for _ in uniprotid])
     # Construct the SQL query
     sql_query2 = f"""
     SELECT *
     FROM "EMBL"
     WHERE uniprotid IN ({placeholders2});
     """

     # Execute the SQL query
     finalDF = pd.read_sql_query(sql_query2, conn, params=uniprotid)
     finalDF = finalDF.merge(uniref50DF, how='outer').drop_duplicates()


     # Close the cursor and connection
     cursor.close()
     conn.close()

     # Return the result as a pandas Series
     return finalDF

def pid2uniref50(pids):
     import pandas as pd
     import sqlite3
     from rotifer.db import ncbi

     i = ncbi.IPGCursor()
     i = i.fetchall(pids)
     i = pd.concat(i)

    # Connect to the SQLite database
     conn = sqlite3.connect('/netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db')
     cursor = conn.cursor()
     pids = i.pid.tolist()

     placeholders = ', '.join(['?' for _ in pids])
     sql_query = f"""
         SELECT *
         FROM EMBL
         WHERE [EMBL-CDS] IN ({placeholders});
     """
     uniref50DF = pd.read_sql_query(sql_query, conn, params=pids)
     uniref50DF['Uniref50'] = 'UniRef50_' + uniref50DF.uniprotid
     # Close the cursor and connection
     cursor.close()
     conn.close()

     # Return the result as a pandas Series
     return uniref50DF




def mview(seqobj,output='sequence.html', background='black', consensus = [100,90,80,70,60], organism=False, find=False):
    '''
    Using Mview to create a HTML file of the sequence algiment. 
    Need to add rule for replace the sequence object description column for organism name when organims is True
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cwd = os.getcwd()
    consensus = ','.join(map(str,consensus))
    f = ''
    if find:
        f = f'-find {find}'
    l2 = '-label2'
    if organism:
        l2 = ''
    if background =='black':
        b = 'black'
        t = 'white'
    else:
        b='white'
        t='black'

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        seqobj.to_file(f'{tmpdirname}/seqfile') 
        Popen(f'mview -in fasta -html head -consensus on -con_threshold {consensus} -label0 {l2} -label4 -label5 -pagecolor {b} -labcolor {t} -textcolor {t} -alncolor {b} -coloring any -con_coloring any {tmpdirname}/seqfile -css off -colormap CLUSTAL -bold -con_colormap CLUSTAL {f} > {cwd}/{output}',
              stdout=PIPE,
              shell=True
              ).communicate()
        
    return f'Consensus file saved on {cwd}/{output}' 

def rpsblast2table(psiblast_output, simple_profiledb=True):
    '''
    Function to read a raw string with a psi blast output and parse it to a table
    ''' 
    import pandas as pd
    import numpy as np

    l = pd.Series(psiblast_output.splitlines())
    f = l.str.startswith(('Query', 'Sbjct', ' Score','>', 'Query='))
    l = l.loc[f].to_frame()
    l.columns = ['text']
    l.loc[l.text.str.startswith("Query="), 'ID'] = l.query('text.str.startswith("Query=")').text.str.split(expand=True)[1]
    l.loc[l.text.str.startswith(">"), 'domain'] = l.query('text.str.startswith(">") == True').text.str.split(expand=True)[0].str.strip('>')
    l.loc[l.text.str.startswith(">PWD"), 'domain'] = l.query('text.str.startswith(">PWD") == True').text.str.split(expand=True)[2].str.strip(',')
    l.loc[l.text.str.startswith(">"), 'source'] = 'allprofiles'
    l.loc[l.text.str.startswith(">PWD"), 'source'] = 'Pfam'
    if simple_profiledb:
        l.loc[l.source=="allprofiles", 'domain'] = l.query('source =="allprofiles"').domain.str.split('.', expand=True)[0].str.split("_", expand=True)[0]
    l.loc[l.text.str.startswith(" Score"), 'evalue'] =  pd.to_numeric(l.query('text.str.startswith(" Score") == True').text.str.split(',',expand=True)[1].str.split('=', expand=True)[1], errors='coerce')
    l.loc[l.text.str.startswith("Query "), 'start'] = pd.to_numeric(l.query('text.str.startswith("Query ") == True').text.str.split(expand=True)[1], errors='coerce')
    l.loc[l.text.str.startswith("Query "), 'end'] = pd.to_numeric(l.query('text.str.startswith("Query ") == True').text.str.rsplit(expand=True, n=1)[1], errors='coerce')
    l.loc[l.text.str.startswith("Sbjct "), 'hit_start'] = pd.to_numeric(l.query('text.str.startswith("Sbjct ") == True').text.str.split(expand=True)[1], errors='coerce')
    l.loc[l.text.str.startswith("Sbjct "), 'hit_end'] = pd.to_numeric(l.query('text.str.startswith("Sbjct ") == True').text.str.rsplit(expand=True, n=1)[1], errors='coerce')
    l.evalue = l.evalue.ffill()
    l.domain = l.domain.ffill()
    l.ID = l.ID.ffill()
    l.source = l.source.ffill()
    l = l.query('text.str.startswith("Query ") or text.str.startswith("Sbjct ")')
    xx = l.groupby(['ID', 'domain', 'evalue', 'source']).agg({'hit_start':'min', 'hit_end':'max', 'start':'min', 'end':'max'}).reset_index()
    return xx

def phobius2table(phobius_output, short=True, add_evalue= 101e-4):
    '''
    Function to read a raw string with a psi blast output and parse it to a table
    Short : Ouputs a short df 
    add_evalue : add a dummy evalue to the short df to later be used in architecture functions.
    ''' 
    import pandas as pd
    l = pd.Series(phobius_output.splitlines())
    l = pd.Series(l)
    f = l.str.startswith(('FT', 'ID'))
    l = l.loc[f].to_frame()
    l.columns = ['text']
    l.loc[l.text.str.startswith("ID"), 'ID'] = l.query('text.str.startswith("ID")').text.str.split(expand=True)[1]
    l.ID = l.ID.ffill()
    l[['phobius', 'start', 'end','prediction']] = l.query('text.str.startswith("FT")').text.str.strip().str.split(expand=True)[[1,2,3,4]]
    l = l.query('text.str.startswith("FT")').iloc[:, 1:6]
    if short:
        l.phobius = l.phobius.replace({"SIGNAL": "SP", "TRANSMEM":"TM"})
        l.prediction = l['prediction'].fillna(l['phobius']).replace({"CYTOPLASMIC.":"IN", "NON":"OUT"})
        l = l.query('prediction in ["TM", "IN","OUT","SP"]')
        l = l[['ID','prediction','start', 'end']].rename({'prediction': 'phobius'}, axis=1)
        if add_evalue:
            l['evalue'] = add_evalue
    return l
def TMprediction(seqobj,
             predictior='phobius',
             cpu=96):
    '''
    Transmembrane predictor it accepts sequence object.
    Only working with phobius, more predictior to be added
    '''

    import tempfile
    from subprocess import Popen, PIPE
    from rotifer.devel.beta.sequence import sequence as sequence
    seqobj = seqobj.copy()

    # Shoul I remove gaps here??
    if predictior == 'phobius':
        with tempfile.TemporaryDirectory() as tmpdirname:
            # temporary save fasta sequence file
            if isinstance (seqobj, sequence):
                seqobj.to_file(f'{tmpdirname}/seqfile') 
                Popen(f'cat {tmpdirname}/seqfile | parallel -j {cpu} --pipe --recstart ">" phobius > {tmpdirname}/out 2> err',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            with open(f'{tmpdirname}/out') as f:
                    phobius_result = f.read()
        
    return (phobius_result) 

def rpsblast(seqobj,
             db=['allprofiles', 'pwld_new_pfam'],
             cpu=96):
    '''
    RPSblast it can accept sequence object. 
    DB can be suplied as python list, values accepted (NIH servers):
    allprofiles
    pwld_new_pfam
    pwld_pfam
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    import os

    if isinstance(db,list):
        dbpath = f'{os.getenv("DATABASES")}/rpsdb/'
        with tempfile.TemporaryDirectory() as tmpdirname:
            for x in db:
                dbx = f'{dbpath}{x}'
                # temporary save fasta sequence file
                if isinstance (seqobj, sequence):
                    seqobj =seqobj.copy()
                    seqobj.to_file(f'{tmpdirname}/seqfile') 
                    Popen(f'cat {tmpdirname}/seqfile | splishrps -d {dbx} -a {cpu} >> {tmpdirname}/out', stdout=PIPE, shell=True).communicate()
                    with open(f'{tmpdirname}/out') as f:
                            blast_r = f.read()
                else:
                    print('please supply a sequence obj')
        return blast_r 
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            # temporary save fasta sequence file
            if isinstance (seqobj, sequence):
                seqobj =seqobj.copy()
                seqobj.to_file(f'{tmpdirname}/seqfile') 
                Popen(f'cat {tmpdirname}/seqfile | splishrps -d {db} -a {cpu} > {tmpdirname}/out', stdout=PIPE, shell=True).communicate()
                with open(f'{tmpdirname}/out') as f:
                        blast_r = f.read()
            else:
                print('please supply a sequence obj')

        return blast_r 

def get_plen(pidlist):
    """
    small function to get plen of proteins, can receive list, pd.Series or a sequence objct
    """
    import pandas as pd
    import rotifer.devel.beta.sequence as rdbs
    if isinstance(pidlist, pd.Series):
        pidlist = pidlist.drop_duplicates().tolist()
        
    if isinstance(pidlist, rdbs.sequence):
        results = pidlist
    else:
        results = rdbs.sequence(pidlist)
    
    results = results.df[['id','length']].rename({'id':'ID', 'length': 'full_length'}, axis=1)
    return results

def domtable(dom_table, seqobj=False):
    import pandas as pd
    import numpy as np
    from rotifer.devel.alpha import gian_func as gf
    copy = dom_table.copy()
    if seqobj:
        dom_table = dom_table.merge(gf.get_plen(seqobj))
    else:    
        dom_table = dom_table.merge(gf.get_plen(dom_table.ID.unique().tolist()))

    dom_table = dom_table[['ID', 'domain', 'start', 'end', 'evalue', 'full_length']].sort_values(['ID', 'start', 'end'])
    dom_table['region'] = np.where(dom_table.ID.shift(+1) == dom_table.ID, dom_table.start - dom_table.end.shift(+1),0)
    dom_table['region2'] =  (dom_table.start - dom_table.region) +1
    dom_table['region3'] =  (dom_table.start - dom_table.region2)
   # unk_central_dom = dom_table.query('region3 >=5')
   # unk_n_term_dom = dom_table.query('region3 ==-1 and region2 >=5')
    dom_table['region4'] = np.where(dom_table.ID.shift(-1) == dom_table.ID,0, dom_table.full_length - dom_table.end)
    dom_table['region5'] = dom_table.start.shift(-1)-1
    dom_table['region5'],dom_table['region6'] = np.where(dom_table.region4 > 0 ,[dom_table.end, dom_table.full_length], [ dom_table.end, dom_table.region5])
    dom_table['region7'] = dom_table.region6 - dom_table.region5
    dom_n = dom_table.drop_duplicates('ID')
    dom_n.end = dom_n.start
    dom_n.start = 1
    dom_n = dom_n[['ID', 'start', 'end']]
    dom_table = dom_table[['ID', 'region5' ,'region6', 'region7']]
    dom_table.columns = ['ID', 'start', 'end', 'rlen']
    dom_n['rlen'] = dom_n.end - dom_n.start
    dom_table = pd.concat([dom_table, dom_n]).sort_values(['ID', 'start', 'end'])
    xx = dom_table.query('rlen >=2')
    xx.domain =['unk']
    dom_file = copy[['ID', 'domain', 'start', 'end']]
    xx['domain'] = 'unk'
    final = pd.concat([dom_file,xx]).sort_values(['ID', 'start','end'])[['ID','domain', 'start','end']]

    return (final)

def draw_architecture(id_list, file, size_median=True, domain_rename=True):

    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import yaml
    import numpy as np
    import seaborn as sns
    import yaml
    import pandas as pd
    from rotifer.devel.beta.sequence import sequence
    import rotifer.devel.alpha.gian_func as gf
    from rotifer.interval import utils as riu
    from rotifer.core import functions as rcf

    svg_dict = {'line_svg' : rcf.findDataFiles(":templates/line.svg"),
                'ANKs' : rcf.findDataFiles(":templates/Ank.svg"),
                'Bps' : rcf.findDataFiles(":templates/Bet_propelers.svg"),
                'HEATs' : rcf.findDataFiles(":templates/HEAT.svg"),
                'LRRs' : rcf.findDataFiles(":templates/LRR.svg"),
                'Pps' : rcf.findDataFiles(":templates/Pps.svg"),
                'Sel1' : rcf.findDataFiles(":templates/Sel1.svg"),
                'TPRs' : rcf.findDataFiles(":templates/TPR.svg")
                }

    #repeats = ['ANKs',  'Bps','HEATs',  'LRRs',  'Pps'  'Sel1',  'TPRs']
    # As reference Nature's standard figure sizes are 89 mm (3.50 inches) wide (single column) and 183 mm (7.20472 inches) wide (double column)
    # If a protein with 5000 aminoacids would fit a whole line of double colum figure, each aa would then use 0.001440944 inches.
    # If I use a PDF with 8 inches and 0.5 margin I would have 7 inches available
    # Therefor we could use this number as conversion factor to draw the domains shapes and rouglth estimate the size to write within it.
    #The font size is typically measured in points. One point is equal to 1/72 of an inch

    # Example : a domain of 300 aa would have a shape of 300 *0.001440944 = 0.4322832 inches
    # If I use font of size 5, each character would 
    # Font size (in inches)= 72/5 ≈0.069 inches
    #For proportional fonts like Arial, the average width of a character is typically about 0.6 times the font size.
    #Average character width (in inches)=0.069×0.6≈0.041 inches

    #Therefore to use calculate the number of characters that fits 300 aa domain:
    # 0.4322832/0.041 == ~10


    # Figure to scale (reccomended to keep 1 to produce final figure in adobe (best scale for publication, but needs to be proper adjusted the names)
    ##### To create PDF, the 4 works preatty good

    seqobj = sequence(id_list)
    seqobj.arch = gf.rpsblast(seqobj)
    tsv =  gf.rpsblast2table(seqobj.arch)
    tsv = tsv.query('evalue <= 0.1')
    seqobj.phobius = gf.TMprediction(seqobj)
    pt = gf.phobius2table(seqobj.phobius).rename({"phobius":"domain"}, axis=1)
    tt = pd.concat([tsv,pt])
    tt.start = tt.start.astype(int)
    tt.end = tt.end.astype(int)
    no = riu.filter_nonoverlapping_regions(tt.query('domain not in ["IN", "OUT"]') ,**{**riu.config['archtsv'], 'maximum_overlap':0.4})
    
    # get the organism name from the fasta description:
    seqobj.df['organism'] = seqobj.df.description.apply(gf.extract_organism_from_description)
    organism_dict = seqobj.df.set_index("id").organism.to_dict()




    scale_figure = 2
    domain_height = scale_figure/10
    font_point = float(0.013837) * 0.6  # 1/72 (size of one point character in inches) * actual size (avg 0.6 the real size)
    aa_scale = float(0.001440944) 
    colors = sns.color_palette("pastel").as_hex() +sns.color_palette("deep").as_hex()  + sns.color_palette('muted').as_hex() + sns.color_palette('colorblind').as_hex()
    shapes = ['box',
              'octagon',
    #          'circle',
              'hexagon',
              'ellipse',
              'diamond']
    #          'doubleoctagon',
    #          'tripleoctagon']


    shapes = pd.DataFrame(shapes, columns=['shapes']).reset_index().rename({'index':'shape_rank'}, axis=1)
    colors = pd.DataFrame(colors, columns=['colors']).reset_index().rename({'index':'color_rank'}, axis=1)
    color_code = shapes.join(colors, how='cross').join(pd.DataFrame(['rounded,filled', 'filled']), how='cross').rename({0:'rounded'}, axis=1).sort_values(['rounded','shape_rank','color_rank'], ascending=[False, True, True]).reset_index(drop=True)
    color_code = color_code.query('~(shapes =="ellipse" and rounded =="rounded,filled")')



    font_size= 4
    d = gf.domtable(no, seqobj=seqobj)
    d = d.rename({"ID":'pid'}, axis=1)

    if domain_rename:
        domain_dict = rcf.loadConfig(rcf.findDataFiles(":data/domain_rename.yaml"))
        domain_dict =  {v:k for k in domain_dict.keys() for v in domain_dict[k] }
        d['domain'] = d['domain'].replace(domain_dict)
    #due Some strange behaivour I will strip the columns domain:
    d['domain'] = d['domain'].str.strip()
    d['size']  = d.end - d.start
    d['esc'] = d['size'] * (aa_scale*scale_figure)
    #d.domain = d.domain.replace(domain_dict['Domain'])
    d['dl'] = d.domain.str.len()
    d['f_space']= round(d.esc/(font_point*font_size))
    d = d[d.f_space >= 1]
    d.f_space = d.f_space.astype(int)
    d['pid_order'] = d.pid.map(d.pid.drop_duplicates().reset_index(drop=True).reset_index().set_index('pid')['index'].to_dict())
    d = d.reset_index(drop=True).reset_index()

    if size_median:
        d.loc[d.domain !='unk', 'esc'] = d.loc[d.domain !='unk'].groupby('domain').esc.transform('median')

    # Adding custoum color and shape for TMs
    ################################################
    ll = d[~d.domain.isin(['unk', 'TM'])]
    ll = ll[~ll.domain.isin(svg_dict.keys())].domain.value_counts().reset_index()
    
    oo = ll.join(color_code)
    ## Adding the TM color and shape
    oo = pd.concat([pd.DataFrame(['TM', 0, 1, 'box', 0, '#FF912B', 'rounded,filled'], index=oo.columns).T, oo])
    ## fillling color and shape for domains not much represented on the ds
    oo.colors = oo.colors.fillna("#808080")
    oo.shapes = oo.shapes.fillna("ellipse")
    oo.rounded = oo.rounded.fillna("filled")
    sd = oo[['index', 'shapes']].set_index('index').shapes.to_dict()
    cd = oo[['index', 'colors']].set_index('index').colors.to_dict()
    rd = oo[['index', 'rounded']].set_index('index').rounded.to_dict()

    ######################################
    #Setting font size based on the shapes
    d['font_size'] = font_size
    d['domain_char'] = d.domain.str.len()
    d['write_size'] = d.domain_char * font_point * font_size
    d['over_write'] = d.esc / d['write_size']
    d.loc[d.over_write < 1, 'font_size'] =   d.loc[d.over_write < 1].font_size * d.loc[d.over_write < 1].over_write
    fsd = d.drop_duplicates('domain').set_index('domain').font_size.to_dict()


    draw_scale=True
    scale_size = 100 * (aa_scale*scale_figure)

    gb = d.groupby('pid_order')
    li =[]
    A = pgv.AGraph()
    #Setting the size of output, based in the number of proteins:
    A.graph_attr['size'] = f'12,{d.pid.nunique()}'
    A.graph_attr['margin'] = '0.5'
    if draw_scale:
        l =['scale_text', 'scale']
        A.add_node('scale_text',
                   label='scale',
                   shape='none',
                   fontsize=font_size +1,
                   fontname='Arial')
        A.add_node('scale',
                   label='100',
                   width=scale_size,
                   fixedsize='true',
                   nodesep=0,
                   fontsize=font_size,
                   #shape='box',
                   image=svg_dict['line_svg'],
                   #style = 'filled',
                   height='0.00001',
                   #fillcolor='black',
                   penwidth=0,
                   labelloc='t',
                   #color='white',
                   labeldistance=1.5)
        li.append(l)

    for x in gb:
        l =[]
        for y,z in x[1].iterrows():
            if len(l) < 1:
                l.append(z['pid'])
                A.add_node(z['pid'],
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="5"></TD></TR>'
                               f'<TR><TD HEIGHT="5"></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="5">{x[1].pid.unique()[0]}</FONT></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial italic" POINT-SIZE="4">{organism_dict[x[1].pid.unique()[0]]}</FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=font_size +1,
                           fontname='Arial')

            if z.domain =='unk':
                A.add_node(z['index'],
                           label='',
                           width=z['esc'],
                           fixedsize='true',
                           nodesep=0,
#                           shape='box',
                           image=svg_dict['line_svg'],
                           height='domain_height',
                           penwidth=0)
                l.append(z['index'])
            elif z.domain in svg_dict.keys():
                A.add_node(z['index'],
                           label=z.domain,
                           width=z['esc'],
                           fixedsize='true',
                           nodesep=0,
                           image=svg_dict[z.domain],
                           fontsize=fsd[z['domain']],
                           height='domain_height',
                           penwidth=0,
                           fontname='Arial')
                l.append(z['index'])
            else:
                A.add_node(z['index'],
                           label=z.domain,
                           style=rd[z['domain']],
                           width=z['esc'],
                           fixedsize='true',
                           height=domain_height,
                           color=cd[z['domain']],
                           fontsize=fsd[z['domain']],
                           fillcolor=cd[z['domain']],
                           shape=sd[z['domain']],
                           nodesep=0,
                           penwidth=1,
                           fontname='Arial')
                l.append(z['index'])
        li.append(l)
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")

    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.02)
    A.draw(file, prog="dot")


def extract_organism_from_description(text):
    import re
    match = re.search(r'\[(.*?)\]', text)
    if match:
        content = match.group(1)
        words = content.split()
        return ' '.join(words[:2])
    return None


def insert_na_rows_at_indexes(df, indexes):
    import pandas as pd
    import numpy as np
    """
    Inserts rows with NaN values into the DataFrame at specified indexes.
    
    Parameters:
    df (pd.DataFrame): The original DataFrame.
    indexes (list of int): A list of indexes at which to insert new rows with NaN values.
    
    Returns:
    pd.DataFrame: A new DataFrame with the NaN rows inserted.
    """
    # Ensure indexes are sorted in descending order to handle multiple insertions
    sorted_indexes = sorted(indexes, reverse=True)

    # Create a DataFrame for the NaN row
    na_row_df = pd.DataFrame([{col: np.nan for col in df.columns}], columns=df.columns)

    # Convert DataFrame to a list of dictionaries for easier manipulation
    df_list = df.to_dict(orient='records')

    for index in sorted_indexes:
        # Insert NaN row into the list at the given index
        df_list.insert(index, na_row_df.iloc[0].to_dict())

    # Convert the list of dictionaries back to DataFrame
    new_df = pd.DataFrame(df_list)

    return new_df




def operon_fig2(df, domain_dict=False,output_file='operon_fig_out.pdf', domain_column='arch', height=0.28, f=2, fontsize=4):
    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import yaml
    import numpy as np
    import pandas as pd
    from rotifer.devel.alpha import gian_func as gf



    def split_domain(df, column=domain_column,fill ='?', domain_dict=domain_dict, strand=False, after=10, before=10, remove_tm=False):
        '''
        It will use the query as anchor to delect the given number of after and before genes to
        create a dataframe of domains.
        One can use a dcitionary of domain to select the ones that would be kept and rename it
        '''
        def reverse_domains(pid):
            def fix_direction(df):
                df = df.copy()
                q = df.query('query ==1').query('strand ==-1').block_id
                df.strand = np.where(df.block_id.isin(q), df.strand * -1, df.strand)
                return df
            pid = fix_direction(pid)
            strand = pid.strand.unique()[0]
            if strand == -1:
                pid.domp = pid.domp.iloc[::-1].to_list()
            return pid.reset_index(drop=True)

        domdf = df.reset_index(drop=True).select_neighbors(strand=strand, after=after, before=before).copy().reset_index(drop=True).query('type == "CDS"')
        domdf['dom'] = domdf[column].str.split('+')
        domdf = domdf.explode('dom')
        domdf['domp'] = pd.Series(np.where(domdf.pid == domdf.pid.shift(1), 1,0))
        domdf.domp = domdf.groupby(['pid','block_id'])['domp'].cumsum()

        if fill:
            domdf.dom = domdf.dom.fillna(fill)

        # removing tm sig
        if remove_tm:
            domdf['tmc'] = np.where(domdf.dom.isin(["TM","SIG"]),0,1)
            only_tm_sig = list(domdf.groupby('pid').tmc.sum().where(lambda x : x==0).dropna().index)
            domdf.dom = np.where(domdf.pid.isin(only_tm_sig), fill, domdf.dom)
            domdf = domdf.query('dom no in ["TM", "SIG"]')

        domdf = domdf.groupby('pid').apply(
                reverse_domains
                ).reset_index(
                        drop=True
                        ).sort_values([
            'nucleotide',
            'block_id',
            'start',
            'end',
            'domp'])

        if not domain_dict:
            return domdf
            
        #domdf = domdf[domdf.dom.isin(list(domain_dict.keys()))]
        domdf.dom = domdf.dom.replace(domain_dict)


        return domdf

    te = df.query('type =="CDS"').copy()
    te = split_domain(te)
    te['shape'] = 'rectangle'
    # Removing duplicates sequencial domain in same protein
    te = te.loc[~((te.pid.shift() == te.pid) & (te.dom.shift() == te.dom))].reset_index(drop=True)

    l = te.drop_duplicates(['block_id','pid'], keep='first').query('strand ==-1').index
    r = te.drop_duplicates(['block_id','pid'], keep='last').query('strand ==1').index
    te.loc[r,'shape'] = 'rarrow'
    te.loc[l,'shape'] = 'larrow'
    te['height'] = np.where(te['shape'].isin(['larrow', 'rarrow']),  height, height/f)
    te['style'] = 'filled' 
    add_rows = te[(te.strand > te.strand.shift(1) )& (te.block_id == te.block_id.shift(1))].index.tolist()
    te = gf.insert_na_rows_at_indexes(te, add_rows)
    te.pid = te.pid.fillna(te.pid.reset_index()['index'])
    te.block_id = te.block_id.fillna(method='ffill')
    te.loc[te.nucleotide.isna(),['dom','domp','shape','style', 'height', 'width']] = ['|', 0, 'point','', 0, 2] 
    te = te.reset_index(drop=True).reset_index()

    gb = te.groupby('block_id')
    li=[]
    A = pgv.AGraph()
   # return gb
    for x in gb:
        l =[]
        for y,z in x[1].iterrows():
            if len(l) < 1:
                l.append(z['block_id'])
                A.add_node(z['block_id'],
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="4">{x[1].block_id.unique()[0]}</FONT></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial italic" POINT-SIZE="4">{x[1].organism.unique()[0]}</FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Arial')

            A.add_node(z['index'],
                       label=z.dom,
                       shape=z['shape'],
                       width=(len(z['dom'])/10) - 0.05,
                       style=z['style'],
                       height=z['height'],
                       color='lightgray',
                       fillcolor='lightgray',
                       fontsize=fontsize)
            l.append(z['index'])

        li.append(l)
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")


    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.02)
    A.draw(output_file, prog="dot")
    return te


