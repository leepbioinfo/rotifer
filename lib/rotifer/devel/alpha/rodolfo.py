#!/usr/bin/env python3
import rotifer
from rotifer import GlobalConfig
logger = rotifer.logging.getLogger(__name__)

def annotate_community(df, column='community', reference='query', ):
    '''
    add or modify communities in a dataframe
    Parameters:
        df: Pandas dataframe
        dataframe to be annotated
        column:
        rows: which rows to annotate
        reference: 
    '''
    # loc to change info
    # or map
    # df.loc[rows, column] = df[reference].map(dict)
    # dict is a result of a df2comm
    # cut communities by cluster size and inputs two tables
    # annotated table and graph
    return

def dataframe_to_community(df, source='query', target='hit', weight='probability', directed=False):
    '''
    Extract communities from a dataframe using the Louvain algorithm.
    '''
    import pandas as pd
    import numpy as np
    import networkx as nx
    import community

    # The following lines are responsible for transforming the input matrix
    # into a symmetric matrix, i.e. an undirected graph
    if not directed:
        df[source], df[target] = np.where(
            df[source] > df[target],
            [df[source], df[target]],
            [df[target], df[source]],
        )
        df = df.sort_values(weight).drop_duplicates([source, target])

    G = nx.from_pandas_edgelist(df.filter([source, target, weight]), source=source, target=target, edge_attr=weight)
    partition = community.best_partition(G, weight=weight)
    d = (
        pd.DataFrame.from_dict(partition, orient="index")
        .reset_index()
        .rename({"index": "nodes", 0: "community"}, axis=1)
    )
    return d

def build_hhsuite_database(
    df,
    sequences,
    group="c80e3",
    members="c80i70",
    output_directory="hhmdb",
    alignment_method='famsa',
    minimum_group_size=3,
    threads=8,
):
    '''
    Build a HH-suite database for all protein groups.

    Parameters:
      df: a Pandas Dataframe
        The members column should list all sequences the user
        wishes to add to the alignment of each group
      sequences: string
        Path to a FASTA file with the sequences under analysis.
      group: string, default 'c80e3'
        Name of the column containing group identifiers
      members: string, default 'c80i70'
        Name of the column storing sequence identifiers
      output_directory: path-list string, default 'hhmdb'
        Name of the directory to save all output files
      minimum_group_size: integer, default 3
        Size of the smallest cluster that must be added
        to the database
      threads: integer, default 10
        Number os CPUs to use while building alignments
    '''
    import os
    import tempfile
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence

    # Create output directory
    if not os.path.exists(f'{output_directory}/aln'):
        os.makedirs(f'{output_directory}/aln')

    # Build MSAs
    for group_cluster in df[group].unique():
        g = df.loc[df[group] == group_cluster,members].drop_duplicates().dropna()
        if g.nunique() >= minimum_group_size:
            aln = sequence(g.tolist(),"accession", name=group_cluster, local_database_path=sequences)
            aln = aln.align(method="famsa")
            aln.to_file(f'{output_directory}/aln/{group_cluster}.a3m','a3m')

    # Build hh-suite database
    with tempfile.TemporaryDirectory() as td:
        cmd = rotifer.GlobalConfig['base'] + f'/share/rotifer/scripts/build_hhdb.sh'
        Popen(f'{cmd} {td}/{output_directory} {output_directory}/aln').communicate()
        for orig in ['cs219','a3m.ordered','hhm.ordered']:
            dest = orig.replace(".ordered","")
            os.rename(f'{td}/{output_directory}.{orig}.ffdata',f'./{output_directory}/{output_directory}_{dest}.ffdata')
            os.rename(f'{td}/{output_directory}.{orig}.ffindex',f'./{output_directory}/{output_directory}_{dest}.ffindex')

def hhblits_all(input_directory, output_directory=None, input_suffix="a3m", output_suffix="hhr"):
    """
    Run all alignments in a directory against a HH-suite database.
    """
    from glob import glob
    if not output_directory:
        output_directory = f'{input_directory}/{output_suffix}'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    for aln in glob(f'{input_directory}/*.{input_suffix}'):
        y = os.path.basename(aln).replace(".{input_suffix}","")
        cmd = f'hhsearch --cpu {cpu} -i {aln} -d {input_directory}/{input_directory} -o {output_directory}/{y}.{output_suffix}'
        Popen(cmd, stdout=PIPE, shell=True).communicate()

def hmmer2aln(hmmout, df, threads, folder):
    



    return df
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
            aravind=('profiledb', count_series)
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


def annotation(seqobj, coordinates, delimiter=True):
    '''
    Receive a sequence object and a list with tuples
    containing [(start, end annotation)], to add
    annotation in alignment
    Example:
        annotation(
        seqobj,
        [(35, 217,
        ' domain 1 (probability 90%, WP_xxxxxx)'),
        (240, 380, 'domain 2')
        ] )
    Function in experimental phase, beware of bugs.
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
    Function to load the result of hmm{search,scan} and obtain
    the results for the full protein.
    If error_lines = False, it skip the results that gave error
    loading to the dataframe.
    Those errors only happens if any field present an string with
    more than two spaces in sequence (mostly found in description field).
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
    import rotifer.interval.utils as riu
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

    hhr_nolr = riu.filter_nonoverlapping_regions(
        tdf.astype({'start': int, 'end': int}).reset_index(),
        end='end',
        start='start',
        reference='query',
        criteria=OrderedDict([('Probab', False), ('region_length', True)]
                             ))
    # ISSUE !!!!!! When the rep seq do not have the first domains
    # it breaks the code as there is no way to map the domain.
    # This below is an improvisation that will overpass it,
    # but with error in annotation!!
    hhr_nolr = hhr_nolr[hhr_nolr.end > hhr_nolr.start]
    return hhr_nolr




def annotate_seqobj(seqobj,
                    df,
                    cnt='profiledb'
                    ):
    '''
    Add annotation to seqobject:
        architecture,
        clusters,
        compact_neighborhood.
    '''
    df = df.copy()
    seqobj = seqobj.copy()
    accs = seqobj.df.id.to_list()
    cn = df[
            df.block_id.isin(df.query('pid in @accs').block_id)
            ].compact_neighborhood(cnt)
    seqobj.df = seqobj.df.merge(df.query('pid in @accs')[
        ['pid',
         'block_id',
         'profiledb',
         'pfam',
         'c80e3',
         'c80i70']
        ].rename(
            {'pid': 'id'}, axis=1
            ).set_index('block_id').join(cn), on='id', how='left')

    return seqobj


def psiblast(acc,
             db='/databases/fadb/nr/nr50',
             cpu=18,
             aln = True,
             num_aln = 500):
    '''
    Psiblast search of a seqobj against nr50 (default),

    Dependencies:
        splishpsi
        blast2table
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

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        if isinstance (acc, sequence):
            acc.to_file(f'{tmpdirname}/seqfile') 
            if aln:
                Popen(f'splishpsi -a {cpu} -num_alignments {num_aln} -in_msa {tmpdirname}/seqfile -d {db} > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            else:
                Popen(f'splishpsi -a {cpu} -num_alignments {num_aln} -i {tmpdirname}/seqfile -d {db} > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            Popen(f'blast2table {tmpdirname}/out > {tmpdirname}/out.tsv',
                  stdout=PIPE,
                  shell=True).communicate()
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
    From a search tabular file, filter by coverage, and evalue to make a aligned seqobj.

    Dependencies:
        splishrps
        rps2arch
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

    #seqobj.df.id = seqobj.df.apply(lambda x: f'{x.id}/{x.hitstart}-{x.hitend}', axis=1)
    return seqobj


def add_arch_to_seqobj(seqobj, db='profiledb', cpu=8):
    '''
    Psiblast it can accept sequence object. 

    Dependencies:
        splishrps
        rps2arch
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


def add_arch_to_df(df):
    '''
    Add architecture and clusters info to a NeighborHooddf
    '''

    import tempfile
    from subprocess import Popen, PIPE
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
        Popen(
                '/home/linuxbrew/git/rotifer/share/rotifer/scripts/scanseqs.sh tmp1 seqfile',
                stdout=PIPE,
                shell=True
                ).communicate()
        info = pd.read_csv(
                'tmp1.c100i100.tsv',
                sep='\t',
                names=['c100i100', 'pid']
                )
        info = info.merge(
                pd.read_csv(
                    'tmp1.c80i70.tsv',
                    sep='\t',
                    names=['c80i70', 'c100i100']
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp1.c80e3.tsv',
                    sep='\t',
                    names=['c80e3', 'c80i70']),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp1.aravind.scan.arch',
                    sep='\t',
                    names=['c100i100', 'profiledb'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp1.pfam.scan.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        df = df.merge(info, how='left')
        os.chdir(cwd)

    return df

def full_annotate(seqobj,
                  progress=True,
                  batch_size=8,
                  mirror=["/am/ftp-genomes","/databases/genomes"],
                  threads=8,
                  after=5,
                  before=5,
                  eukaryotes=False):
    '''
    Annotate a seqobj with its own genome neighborhood.
    '''
    from rotifer.db import ncbi 
    from rotifer.devel.alpha.rodolfo import add_arch_to_df
    gnc = ncbi.GeneNeighborhoodCursor(
            progress=progress,
            batch_size=batch_size,
            mirror=mirror,
            threads=threads,
            after=after,
            before=before,
            eukaryotes=eukaryotes)
    seqobj.ndf = add_arch_to_df(gnc.fetchall(seqobj.df.id.to_list()))
    #seqobj.ndf = seqobj.ndf.drop_duplicates(['start', 'end', 'nucleotide'])
    return seqobj

def padding_df(df):
    cdf = df.copy()
    c = df.columns
    pad_col_name=[]
    for x in c:
        cdf[x] = cdf[x].fillna('').astype(str)
        w = cdf[x].str.len().max()
        cdf[x] = cdf[x].str.pad(width =w)
        pad_col_name.append(x.rjust(w))
    cdf.columns = pad_col_name    
    return cdf


def alnxaln(seqobj, clustercol = 'c50i0', minseq=10):
    import os
    import tempfile
    import pandas as pd
    from rotifer.io import hhsuite

    '''
    Compares a seqobj within its clusters using hhalign.

    Dependencies:
        hhalign
        bash
    '''
    
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

def split_by_model(seqobj=None, df=False, evalue=1e-3, arch=None, coverage=50, method='hmmer', ldbp=None):
    '''
    Using a hmm model to split your sequence object to match only the model region.
    If more than one match in one protein, it will split the match in two sequence.
    From a search tabular file, filter by coverage, and evalue to make a aligned seqobj.

    Dependencies:
        hmmersearch
        hmmer2table
        domain2architecture
        architecture2table
        splishrps
        rps2arch
    '''

    import os
    import pandas as pd
    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    
    if coverage > 1:
        coverage = coverage/100

    #seqobj =seqobj.copy()

    if method == 'hmmer':

        df = df.query(f'evalue <= {evalue} and cov > {coverage}')
       # with tempfile.TemporaryDirectory() as tmpdirname:
           # seqobj.to_file(f'{tmpdirname}/seqfile') 
           # Popen(f'hmmsearch {model} {tmpdirname}/seqfile > {tmpdirname}/hmmsearch_r',
           #       stdout=PIPE,
           #       shell=True).communicate()
           # Popen(f'hmmer2table {tmpdirname}/hmmsearch_r |domain2architecture |architecture2table > {tmpdirname}/hmmsearch_table',
           #       stdout=PIPE,
           #       shell=True).communicate()
           # t = pd.read_csv(f'{tmpdirname}/hmmsearch_table', sep="\t")
        seqobj = sequence(df.sequence.drop_duplicates().to_list(), local_database_path=ldbp)
        seqobj.df = seqobj.df.merge(df.rename({'sequence':'id'},axis=1), how='left')
        seqobj.df.end = seqobj.df.eend.fillna(seqobj.df.length)
        seqobj.df.start = seqobj.df.estart.fillna(1)
        seqobj.df.sequence = seqobj.df.apply(lambda x: x.sequence[int(x.estart):int(x.eend)], axis=1)
        seqobj.df['id'] = seqobj.df['id'] + '/' + seqobj.df['estart'].astype(str) + '-' + seqobj.df['eend'].astype(str)

    else:

        if df:
            df = df.query(f'evalue <= {evalue} and querycoverage > {coverage}')
            seqobj = seqobj.filter(keep=df.hit.drop_duplicates().to_list())
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

    return (seqobj)

def remove_redundancy(seqobj, identity =80, coverage = 70):
    from rotifer.devel.beta.sequence import sequence as sequence
    import pandas as pd
    s = seqobj.copy()
    s.add_cluster(coverage=coverage, identity=identity, inplace=True)
    s = s.filter(keep=s.df[f'c{coverage}i{identity}'].drop_duplicates().to_list())
    return(s)


