#!/usr/bin/env python3

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0]))), "lib"))
import rotifer.core.cli as corecli
import rotifer.core.functions as rcf
from rotifer.core.log import log
from rotifer.neighborhood import neighborhood
from rotifer.alchemy.connect import clickhouse
from rotifer.core import GlobalConfig
from rotifer.db.ncbi import NcbiConfig

import subprocess
import pandas as pd
import logging
from multiprocessing import Process, Manager
import datetime
from Bio import SeqIO
from Bio import Entrez
import time
import shutil
import numpy as np
import warnings
from os.path import expanduser
from datetime import datetime as dt
import argcomplete
from clickhouse_driver import Client
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
import random, string
from clickhouse_sqlalchemy import Table, make_session, get_declarative_base, types, engines
#from pandas.compat import StringIO
from io import StringIO

from tqdm import tqdm

warnings.filterwarnings('ignore')

__version__ = 0.62
__authors__ = 'Gilberto Kaihami, Aureliano Guedes, Gabriel Hueck, Robson F. de Souza'
__rdoc__ = '''
DESCRIPTION:
Get gene neighborhood
'''

########
# New features
########
# In version 0.60
# Added rdoc

# In version 0.59
# Added minimal intergenic distance (mid) option

# In version 0.58
# Removed some redundant code
# Fixed duplicated output

# In version 0.56
# Added multiple epost

# In version 0.55
# Modified genomes clickhouse added timestamp

# In version 0.54
# Fixed circular genome problem

# In version 0.53
# Debug dataframe option

# In version 0.52
# Removed redundant code when collecting ipgs

# In version 0.50
# Fixed minor error with the query
# Added a progress bar

# In version 0.47
# Changed how to retrieve neighborhood => 50-100x faster

# In version 0.43:
# Fixed the genbank problem not checking the DB
# Started the docstring

# Fixed ncRNA query and other regulatory RNA regions

# In the version 0.40:
# Using clickhouse to insert the gbff in the database
# Fast feature2df

# TODO: Still need to put efetch/esearch feature

def parse_cli():

    parser = corecli.parser(description = 'Get gene neighborhood')

    parser.add(rconfig = ':cli.acc')

    parser.add(long_arg = '--api_key',
               dest = 'api_key',
               helper = 'Set NCBI API key',
               default = '')

    parser.add(long_arg = '--email',
               dest = 'email',
               helper = 'Set email for NCBI API',
               default = '')

    parser.add(long_arg = '--search_method',
               short_arg = '-sm',
               dest = 'search',
               helper = 'Select a search method, local or www',
               default = 'www')

    parser.add( short_arg = '-a',
                long_arg ='--above',
                dest = 'above',
                helper = 'Rows above: maximum number of neighbors upstream of target loci (Default: 3)',
                default = 3)

    parser.add( short_arg = '-b',
                long_arg = '--below',
                dest = 'below',
                helper = 'Rows below: maximum number of neighbors downstream of target loci (Default: 3)',
                default = 3)

    parser.add( long_arg = '--progress',
                short_arg = '-p',
                dest = 'progress',
                helper = 'Show progress bar',
                action = 'store_true'
                )

    parser.add( long_arg = '--minimal_intergenic_distance',
                short_arg = '-mid',
                dest = 'distance',
                helper ='''Minimal intergenic distance between two queries. (Default: 0)
Example:
if the distance is set 0 (default) and -a and -b options are set 0, this setting can be thought as continous segments
This will generate one block,
block 1:
--> query_A
--> query_B
if the distance is set -1 (at least one gene overlapping between each block),
The same query will generate two blocks
block 1:
--> query_A
block 2:
--> query_B
if the distance is set with positive values (higher than 0) it means the minimun distance between two blocks.
Like, -mid 2, and -a 0, -b 0.
block 1:
--> query_A
--  protein_1
--  protein_2
--> query_B
               ''',
                default = 0,
                )

    parser.add(long_arg = '--missing',
               dest = 'missing',
               helper = 'Basename to save list of protein accessions with no IPG report.',
               short_arg = '-m',
               default = None,
               arg_type = str)

    # Now we have the writer columns can be added and removed

    parser.add(short_arg = '-v', long_arg = '--verbose',
               dest = 'verbose',
                       action = 'count')

    parser.add(short_arg = '-t',long_arg = '--threads',
               dest = 'threads',
                        helper = 'Number of threads (Default: 3)',
                        default = 3)

    parser.add(short_arg = '-of',
               long_arg = '--outformat',
               dest = 'outformat',
                        helper = 'Output format (table/gi2operon) [Default: table]',
                        default = 'table')

    parser.add(short_arg = '-y',
               long_arg ='--header',
               action='append',
               dest = 'addinfo',
               helper = 'Add more information to gi2operon header. The additional information must be a valid collumn in the table format\n\
(example: -y assembly, -y classification -y seq_type)',
               default = [])

    parser.add( long_arg = '--hit_method',
                short_arg = '-hm',
                dest = 'hit',
                default = 'best',
                arg_type = str,
                helper = 'Select hit method (default: best) [best/any]. "Best" get the best hit, any will get any hit present in the db',
                action = "store"
                )

    parser.add( long_arg = '--assembly',
                short_arg = '-asm',
                dest = 'assembly',
                default = [],
                arg_type = str,
                helper = 'Filter by assembly',
                action = 'append'
                )

    parser.add(long_arg = '--ipg',
               short_arg = '-ipg',
               dest = 'ipg',
               default = None,
               arg_type = str,
               helper = 'Basename for tables with IPG report data downloaded from the NCBI.'
               )

    argcomplete.autocomplete(parser)

    args = parser.parse_args()
    return args

def verbose_msg(message = ''):
    now = dt.now().strftime('[%D %H:%M:%S]')
    if isinstance(message, list):
        message = ' '.join(message)
    sys.stderr.write('## {0} {1}\n'.format(now, message))
    sys.stderr.flush()

def nuc2gbk(gbks, verbose = False, api_key = ''):
    '''
    Get genbank file from nucleotide accession using Entrez
    ----------
    PARAMETERS:
    gbks:    list of nucleotides
    verbose: verbose level
    api_key: Entrez API Key
    '''
    for nucleotide in gbks:
        if verbose:
            verbose_msg('Downloading {0} using efetch'.format(nucleotide))
        try:
            to_save = Entrez.efetch(db = 'nuccore',
                                    rettype = 'gbwithparts',
                                    retmode = 'text',
                                    id = nucleotide,
                                    api_key = api_key)

            out_handle = open(GlobalConfig['cache'] + '/' + nucleotide + '.gbff', "w")
            out_handle.write(to_save.read())
            out_handle.close()
            to_save.close()
        except:
            if verbose:
                if isinstance(nucleotide, str):
                    sys.stderr.write(nucleotide+'\n')

def asm2ftp(gbks, verbose = False,progress = False, position = 0, threads = 0):
    '''
    Get genbank file from a list of asssemblies using ftp
    ----------
    PARAMETERS:
    gbks:    list of assemblies
    verbose: verbose level
    '''

    # Progress messages
    if verbose:
        ngbks = len(gbks)
        verbose_msg(f'asm2ftp: Batch {threads}, {position-1} threads, downloading {ngbks} assembly GBFF files...')
    if progress:
        gbks = tqdm(gbks, position = position, leave=False, desc = f'Downloading genomes thread {threads}')

    # Download assemblies
    burl = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/'
    for ele in gbks:
        got_gbks_folder = [x.replace('.gbff', '').replace('.tsv','') for x in os.listdir(GlobalConfig['cache'])]
        if ele not in got_gbks_folder:
            if verbose:
                verbose_msg('Downloading {0} using ftp'.format(str(ele)))
            tries_download = 0

            while tries_download < 15:
                try:
                    #gca, f1, f2, f3 = ele[0:3], ele[4:7], ele[7:10], ele[10:13]
                    ftppath = "/".join([ burl, ele[0:3], ele[4:7], ele[7:10], ele[10:13] ])
                    cmd = 'curl -l ' + ftppath + '/* 2> /dev/null | grep ' + ele + ' 2> /dev/null'
                    sub2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
                    genome = sub2.communicate()[0].decode('utf-8').replace('\n', '')
                    cmd = 'curl -l ' + ftppath + '/* 2> /dev/null | grep ' + ele + ' | parallel curl -l ' + ftppath + '/{}/ 2> /dev/null'
                    sub1 = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
                    features = sub1.communicate()[0].decode('utf-8').splitlines()
                    for feature in features:
                        file2save = GlobalConfig['cache'] + '/' + ele
                        if file2save+'.gbff' not in os.listdir(GlobalConfig['cache']):
                            if '_genomic.gbff' in feature and '_cds' not in feature and '_rna_' not in feature:
                                download = 'curl ' + '/'.join([ftppath, genome, feature]) + ' --output ' + file2save + '.gbff.gz 2> /dev/null'
                                subprocess.Popen([download], shell = True).wait()
                                subprocess.Popen(['gunzip -f {0}.gbff.gz'.format(file2save)], shell = True).wait()
                    break
                except:
                    tries_download +=1
                    time.sleep(0.3)

def long_routine(queries, df_accs, df_accs2, verbose, progress, position, threads, **kwargs):
    ## Try using epost if fail except using efetch
    if verbose:
        verbose_msg(f'Downloading IPG reports! Batch {position}, {threads} threads, {len(acc)} accessions')

    # Prepare rotifer.db.ncbi.ncbi object
    from rotifer.db.ncbi import ncbi as ncbiClass
    from rotifer.core import GlobalConfig

    # Get ipgs for the current batch of 200 queries
    ncbi = ncbiClass()
    ncbi.submit(queries)
    try:
        ipgs = ncbi.read('ipg')
    except:
        sys.stderr.write(f'There was an error downloading IPG reports: {sys.exc_info()[1]}')
        sys.exit()

    # Mark missing data as 'N/A' and choose the last query in each IPG as its representative
    missing_accessions = ncbi.missing()
    ipgs = ipgs.fillna('N/A')
    ipgs['assembly_modified'] = ipgs.apply(lambda x: x.assembly if x.assembly != 'N/A' else x.nucleotide, 1)
    if isinstance(ipgs['representative'].values[0], np.float64):
        sys.stderr.write('No results\n')
        sys.exit()

    # Best genomes from last hit?
    first_ocurrence = ipgs[(ipgs['query'] == 1) & (ipgs['nucleotide'] != 'N/A')]
    first_ocurrence = first_ocurrence.drop_duplicates(['representative'], keep='first')
    second_case = ipgs[(ipgs['representative'].isin(first_ocurrence['representative'].values)) & (ipgs['assembly_modified'].isin(first_ocurrence['assembly_modified'].values))]

    # Hit2
    r2 = ipgs[~( ipgs['representative'].isin(first_ocurrence['accession'].values) )]
    r2 = r2.replace('N/A', 'No nucleotide')
    r2 = r2[r2['assembly_modified'] != 'No nucleotide']

    try:
        r2 = r2.drop_duplicates(['id'], keep = 'first')
    except:
        pass

    # This df contains Sequences with the same acc and different acc
    df2get = pd.concat([second_case, r2])
    df2get = df2get.replace('', 'No asm')
    df2get['assembly'] = df2get['assembly'].replace('N/A', 'No asm')

    # Drop nucleotides with nothing
    df2get = df2get[df2get['nucleotide'] != 'N/A']

    # Connect to clickhouse database to check for cached data
    if verbose:
        nipgs    = len(ipgs)
        nmissing = len(missing)
        nfound = len(found)
        verbose_msg(f'Batch {position}, {nipgs} IPG rows fetched, {nfound} accessions found, {nmissing} missing: connecting to clickhouse...')

    # GBKs in the DB
    uri = 'clickhouse://default:@localhost/rotifer'
    engine = create_engine(uri)
    session = make_session(engine)
    metadata = MetaData(bind=engine)
    metadata.reflect(bind = engine)
    gbksInLocalDB = rcf._flatten([list(x) for x in engine.execute('''SELECT DISTINCT nuc_asm from genomes''').fetchall()])

    #######
    df2get_asms = list(df2get[df2get['assembly'] != 'No asm']['assembly'].unique())
    df2get_nucs = list(df2get[df2get['assembly'] == 'No asm']['nucleotide'].unique())
    missingAssemblies = [x for x in df2get_asms if x not in gbksInLocalDB]

    #######
    if verbose:
        nipgs    = 0 if ipgs.empty else len(ipgs)
        verbose_msg(f'Batch {position}, {nipgs} IPG rows fetched, {len(found)} accessions found, {len(missing)} missing, downloading {len(missingAssemblies)} GBFF files...')

    # Download genomes to cache
    #asm2ftp(missingAssemblies, gbks_folder, verbose, progress = progress, position = position+threads+1, threads = position)
    ncbi.reset() # Erase previous data
    ncbi.submit(missingAssemblies)
    gbksDownloaded = ncbi.read('genomes', concat=True, parser=False, tempfile=False, delete=False).files

    # Missing assemblies
    #got_gbks_folder = [ x.replace('.gbff.gz', '') for x in gbksDownloaded ]
    #missing_asms = [x for x in df2get_asms if x not in got_gbks_folder]
    #missing_asms = [x for x in missing_asms if x not in got_gbks]
    missing_asms = ncbi.missing()

    find_nuc = [x for x in df2get_nucs if x not in got_gbks_folder]
    find_nuc = [x for x in find_nuc if x not in got_gbks]

    try:
        find_nuc = find_nuc + missing_asms
    except:
        if verbose:
            verbose_msg('Error find_nuc')
            verbose_msg(f'find_nuc {find_nuc}')
            verbose_msg(f'missing_asms, {missing_asms}')

    if find_nuc:
        nuc2gbk(find_nuc, gbks_folder, api_key = api_key)

    if df_accs.empty:
        df_accs = df2get
    else:
        df_accs = pd.concat([df_accs, df2get])

    # Insert into Database raw table
    df_accs_all = ipgs[['id', 'accession', 'assembly', 'nucleotide', 'representative', 'assembly_modified']] # copy for all accs
    df_accs_all = df_accs_all[df_accs_all['assembly_modified'] != 'N/A']
    if df_accs2.empty:
        df_accs2 = df_accs_all[['acc', 'representative', 'nucleotide']]
        df_accs2['representative'] = df_accs2['representative'].astype('category')
    else:
        _ = df_accs_all[['acc', 'representative', 'nucleotide']]
        _['representative'] = _['representative'].astype('category')
        df_accs2 = pd.concat([df_accs2, _])

    return (df_accs, df_accs2)

def _multiple_search(subs, verbose, ls, position, threads, progress = False):
    tot = len(subs)

    # Loop over sets of accessions (subs)
    df_accs = pd.DataFrame()
    df_accs2 = pd.DataFrame()
    w = 0
    try:
        # Get ipgs
        if progress:
            for sub in tqdm(subs, position = position, leave=False, desc = f'Search ipg thread {position}'):
                # long routine here
                df_accs,df_accs2 = long_routine(sub, df_accs, df_accs2, verbose,
                                   progress = progress, position = position, threads = threads)
        else:
            for sub in subs:
                if verbose:
                    verbose_msg(f'There are {tot-w} epost left, total is {tot}, thread {position}')
                    w += 1
                df_accs,df_accs2 = long_routine(sub, df_accs, df_accs2, verbose,
                                   progress = progress, position = position, threads = threads)
    except:
        pass

    if df_accs.empty:
        pass
    else:
        ls.append((df_accs,df_accs2))


def acc2asm2gbk(accs, verbose, df_accs, threads = 3, api_key = '', progress = False):
    '''
    From a acc list get genbank
    -----------
    PARAMETERS:
    accs:
    verbose:
    df_accs:
    threads:
    api_key:
    '''
    # Split acc in a list of size s

    #### TODO remove useless stuf here

    s = 400
    accs2find = [acc for acc in accs if acc not in df_accs['acc'].values] # Ok
    subs = [accs2find[x:x+s] for x in range(0, len(accs2find), s)]        # Ok

    size = len(subs)//threads if len(subs)//threads > 0 else 1
    sub_list = [subs[x:x+size] for x in range(0, len(subs), size)]

    manager = Manager()
    epost_list = manager.list()

    jobs = []
    for x in range(0,len(sub_list)):
        p = Process(target = _multiple_search, args = (sub_list[x], verbose, epost_list, x, threads, progress))
        p.start()
        jobs.append(p)

    try:
        for p in jobs:
            p.join()
    except KeyboardInterrupt:
        for p in jobs:
            p.terminate()
        sys.exit(2)

    if epost_list:
        df_accs = pd.concat([x for x,y in epost_list])
        df_accs2 = pd.concat([y for x,y in epost_list])
    else:
        df_accs = pd.DataFrame()
        df_accs2 = pd.DataFrame()
    verbose_msg(f'acc2asm2gbk: {threads} threads, accessions: {accs}')

    # sys.stderr.write(df_accs.to_string())

    return (df_accs, df_accs2)

def _multiple_find_acc(sub_fi, gbk, send_end):
    results = []
    for f in sub_fi:
        df = pd.read_csv(os.path.join(GlobalConfig['cache'],'acc/'+f), sep = '\t')
        df = df[df['assembly'] == gbk]['acc'].values
        results.extend(df)
    send_end.send(results)

def features2df(gbk, pair_acc = '', verbose = False):
    '''
    Parameters:
    gbk:  assembly or nucleotide
    df:   empty dataframe
    pair_acc
    '''

    # Create engine
    uri = 'clickhouse://default:@localhost/rotifer'
    engine = create_engine(uri)
    session = make_session(engine)
    metadata = MetaData(bind=engine)
    metadata.reflect(bind = engine)
    date = datetime.date(int(datetime.datetime.today().strftime('%Y')), int(datetime.datetime.today().strftime('%m')), int(datetime.datetime.today().strftime('%d')))

    # Check if gbk is present in the database
    if not engine.execute(f"select DISTINCT nuc_asm FROM genomes where assembly = '{gbk}'").fetchone():
        if verbose:
            verbose_msg(f'## loading {gbk}')

        parser_start_time = time.time()
        fts = SeqIO.index(os.path.join(GlobalConfig['cache'],'gbff')+'/'+gbk+'.gbff', "genbank")
        parser_end_time = time.time()
        if verbose:
            verbose_msg(f'## total time to load: {parser_end_time-parser_start_time}')
        not_types = []

        tot = 0
        for k,rec in fts.items():
            seq_type = 'chromosome'
            gbk_id = k
            topology = rec.annotations['topology']
            organism = rec.annotations['organism']
            taxonomy = '; '.join(rec.annotations['taxonomy'])
            try:
                assembly = [x.split(':')[-1] for x in rec.dbxrefs if 'Assembly' in x][0]
            except:
                assembly = '.'
            if not assembly:
                assembly = '.'

            if rec.features:
                res = []
                ft_order = {}

                for (idx, ft) in enumerate(rec.features):
                    if 'plasmid' in ft.qualifiers:
                        seq_type = 'plasmid'
                    if ft.type not in not_types:
                        feature_type = ft.type
                        sub_fts = ft.qualifiers
                        pid = sub_fts['protein_id'][0] if 'protein_id' in sub_fts.keys() else ''
                        if pid == '' and 'CDS' in ft.type:
                            feature_type = 'PSE'

                        # this part is new
                        if feature_type in ft_order.keys():
                            ft_order[feature_type] +=1
                        if feature_type not in ft_order.keys():
                            ft_order[feature_type] = 1

                        locus_tag = sub_fts['locus_tag'][0] if 'locus_tag' in sub_fts.keys() else '.'
                        location = ft.location
                        start = str(location.start+1).replace('>', '').replace('<', '')
                        end = str(location.end).replace('>', '').replace('<', '')
                        try:
                            strand = str(location.strand)
                        except:
                            strand = '.'

                        gene = sub_fts['gene'][0] if 'gene' in sub_fts.keys() else '.'
                        plen = len(ft.qualifiers['translation'][0]) if 'translation' in ft.qualifiers.keys() else 0
                        product = sub_fts['product'][0] if 'product' in sub_fts.keys() else ''
                        res.append({'internal_id':idx, 'nucleotide': gbk_id,
                                    'assembly': assembly,
                                    'topology': topology, 'start': str(start),
                                    'end':str(end), 'strand': strand,
                                    'pid': pid, 'type': feature_type,
                                    'plen': str(plen), 'locus': locus_tag,
                                    'seq_type': seq_type, 'gene':gene,
                                    'product': product, 'organism': organism,
                                    'taxonomy': taxonomy,
                                    'feature_order':ft_order[feature_type],
                                    'nuc_asm': gbk,
                                    'date': date})

                # write res here
                client = Client('localhost')

                # Insert
                insert_start_time = time.time()
                if len(res) > 0:
                    if verbose:
                        verbose_msg(f'## Inserting into clickhouse DB {gbk}')
                    client.execute('Insert into rotifer.genomes values', res)
                    insert_end_time = time.time()
                    if verbose:
                        verbose_msg(f'## Total time to insert {insert_end_time-insert_start_time}')

def _multi_ft2df(ls_gbks = [], verbose = False):
    for gbk in set(ls_gbks):
        if verbose:
            verbose_msg('Ft2df of {0}'.format(gbk))
        try:
            df = pd.DataFrame()
            features2df(gbk, df, verbose = verbose)
        except:
            if verbose:
                _ = 'Error ft2df {0}'.format(gbk)
                verbose_msg(str(_))

###################################

def block_intervals(accs = '', above = 3, below = 3,
              distance = 0, block_id = 0, of = 'table',
              conn = '',
              nucleotide = '',
              nuc_asm = '',
              other_info = [],
              original_acc = []
              ):

    acc_formated = ','.join(["'" + x +"'" for x in accs])

    accs_indexes = conn.execute(f"""SELECT feature_order from genomes
                                where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                pid in ({acc_formated}) ORDER BY feature_order""").fetchall()
    accs_indexes= sorted(list(zip(*accs_indexes))[0])

    max_index, min_index = conn.execute(f"""SELECT
                                        max(feature_order), min(feature_order) from genomes
                                where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                type = 'CDS'""").fetchone()

    max_internal_id, min_internal_id =  conn.execute(f"""SELECT max(internal_id), min(internal_id)
                                         from genomes
                                where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}'""").fetchone()

    topology = conn.execute(f"""SELECT distinct(topology) from genomes
                                where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                type = 'CDS'""").fetchone()[0]
    intervals = []
    for value in accs_indexes:

        find_up = value - above
        find_down = value + below
        if intervals:
            if intervals[-1][1] - find_up+1 >= distance:
                if above+below+1 >= max_index:
                    intervals[len(intervals)-1] = [0, max_index]
                else:
                    intervals[len(intervals) -1 ] = [intervals[len(intervals)-1][0], find_down]
            else:
                if above + below +1 >= max_index:
                    intervals.append([0, max_index])
                else:
                    intervals.append([find_up, find_down])
        else:
            if above +below +1 >= max_index:
                intervals.append([0, max_index])
            else:
                intervals.append([find_up, find_down])
    for up_index, down_index in intervals:
        sub_df = pd.DataFrame()
        get_more_up = ''
        get_more_down = ''
        if up_index < min_index:
            if topology != 'linear':
                get_more_up = max_index+ up_index
            up_index = min_index

        if down_index > max_index:
            if topology != 'linear':
                get_more_down = down_index - max_index
            down_index = max_index

        if get_more_up:
            _ = conn.execute(f"""
                                SELECT internal_id from genomes
                                WHERE nuc_asm = '{nuc_asm}' and
                                nucleotide = '{nucleotide}' and
                                type = 'CDS'
                                and feature_order in ({max_index}, {get_more_up})
                                and type in ('CDS', 'PSE', 'tRNA', 'ncRNA', 'rRNA', 'tmRNA')
                                  """).fetchall()
            _ = [x[0] for x in _]
            if len (_) == 1:
                max_order = min_order = _[0]
                min_order = max_internal_id
            else:
                max_order, min_order = sorted(_)
                min_order = max_internal_id

            q = (conn.execute(f"""SELECT
                                          internal_id,nucleotide,
                                                   nuc_asm, topology, start, end,
                                                   strand,pid, type, plen, locus,
                                                   seq_type, gene, product,organism,
                                                   taxonomy
                                  FROM genomes where nuc_asm = '{nuc_asm}'
                                                      and nucleotide = '{nucleotide}'
                                                      and internal_id >= {max_order}
                                                      and internal_id <= {min_order}
                                  and type in ('CDS', 'PSE', 'tRNA', 'ncRNA', 'rRNA', 'tmRNA')
                                  ORDER BY internal_id""").fetchall())
            sub_df = pd.DataFrame(q,
                                        columns = ['internal_id','nucleotide',
                                                   'assembly', 'topology', 'start', 'end',
                                                   'strand','pid', 'type', 'plen', 'locus',
                                                   'seq_type', 'gene', 'product','organism',
                                                   'classification', ])

        if not sub_df.empty:
            _ = conn.execute(f"""
                                SELECT internal_id from genomes
                                WHERE nuc_asm = '{nuc_asm}' and
                                nucleotide = '{nucleotide}' and
                                type = 'CDS'
                                and feature_order in ({up_index}, {down_index})
                                  """).fetchall()
            _ = [x[0] for x in _]
            if len (_) == 1:
                max_order = min_order = _[0]
            else:
                max_order, min_order = sorted(_)

            if get_more_up:
                max_order = min_internal_id

            if get_more_down:
                min_order = max_internal_id
            q = (conn.execute(f"""SELECT
                                          internal_id,nucleotide,
                                                   nuc_asm, topology, start, end,
                                                   strand,pid, type, plen, locus,
                                                   seq_type, gene, product,organism,
                                                   taxonomy
                                  FROM genomes where nuc_asm = '{nuc_asm}'
                                                      and nucleotide = '{nucleotide}'
                                                      and internal_id >= {max_order}
                                                      and internal_id <= {min_order}
                                  AND type in ('CDS', 'PSE', 'tRNA', 'ncRNA', 'rRNA', 'tmRNA')
                                  ORDER BY internal_id""").fetchall())
            t = pd.DataFrame(q,
                                        columns = ['internal_id','nucleotide',
                                                   'assembly', 'topology', 'start', 'end',
                                                   'strand','pid', 'type', 'plen', 'locus',
                                                   'seq_type', 'gene', 'product','organism',
                                                   'classification', ])
            sub_df = pd.concat([sub_df, t])

        else:
            _ = conn.execute(f"""
                                SELECT internal_id from genomes
                                WHERE nuc_asm = '{nuc_asm}' and
                                nucleotide = '{nucleotide}' and
                                type = 'CDS'
                                and feature_order in ({up_index}, {down_index})
                                  """).fetchall()

            _ = [x[0] for x in _]

            if len (_) == 1:
                max_order = min_order = _[0]
            else:
                max_order, min_order = sorted(_)

            if get_more_up:
                max_order = min_internal_id

            if get_more_down:
                min_order = max_internal_id

            q = (conn.execute(f"""SELECT
                                          internal_id,nucleotide,
                                                   nuc_asm, topology, start, end,
                                                   strand,pid, type, plen, locus,
                                                   seq_type, gene, product,organism,
                                                   taxonomy
                                  FROM genomes where nuc_asm = '{nuc_asm}'
                                                      and nucleotide = '{nucleotide}'
                                                      and internal_id >= {max_order}
                                                      and internal_id <= {min_order}
                                  AND type in ('CDS', 'PSE', 'tRNA', 'ncRNA', 'rRNA', 'tmRNA')
                                  ORDER BY internal_id""").fetchall())

            sub_df = pd.DataFrame(q,
                                        columns = ['internal_id','nucleotide',
                                                   'assembly', 'topology', 'start', 'end',
                                                   'strand','pid', 'type', 'plen', 'locus',
                                                   'seq_type', 'gene', 'product','organism',
                                                   'classification'])

        if get_more_down:

            _ = conn.execute(f"""
                                SELECT internal_id from genomes
                                WHERE nuc_asm = '{nuc_asm}' and
                                nucleotide = '{nucleotide}' and
                                type = 'CDS'
                                and feature_order in ({min_index},{get_more_down})
                                  """).fetchall()

            _ = [x[0] for x in _]
            if len (_) == 1:
                max_order = min_order = _[0]
                max_order = min_internal_id
            else:
                max_order, min_order = sorted(_)
                max_order = min_internal_id

            q = (conn.execute(f"""SELECT
                                          internal_id,nucleotide,
                                                   nuc_asm, topology, start, end,
                                                   strand,pid, type, plen, locus,
                                                   seq_type, gene, product,organism,
                                                   taxonomy
                                  FROM genomes where nuc_asm = '{nuc_asm}'
                                                      and nucleotide = '{nucleotide}'
                                                      and internal_id >= {max_order}
                                                      and internal_id <= {min_order}
                                  AND type in ('CDS', 'PSE', 'tRNA', 'ncRNA', 'rRNA', 'tmRNA')
                                  ORDER BY internal_id""").fetchall())

            if not sub_df.empty:
                t = pd.DataFrame(q,
                                            columns = ['internal_id','nucleotide',
                                                       'assembly', 'topology', 'start', 'end',
                                                       'strand','pid', 'type', 'plen', 'locus',
                                                       'seq_type', 'gene', 'product','organism',
                                                       'classification'])

                sub_df = pd.concat([sub_df, t])


        if not sub_df.empty:
            block_id +=1
            sub_df['block_id'] = sub_df.shape[0]*[block_id]

            sub_df['query'] = sub_df['pid'].map(lambda x: 1 if x in accs else 0)
            if block_id == 1:
                print_header = True
            else:
                print_header = False
            header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene modified product organism classification'.split()
            try:
                dc = {k:v for k,v in zip(original_acc['acc'].values, original_acc['representative'].values)}
                sub_df['modified'] = sub_df['pid'].map(lambda x: dc[x] if x in dc.keys() else '.')
                sub_df['modified'] = sub_df.apply(lambda x:'.' if x['pid'] == x['modified'] else x['modified'],1 )
            except:
                sub_df['modified'] = '.'

            sub_df = sub_df.drop_duplicates()
            neighborhood.writer(sub_df[header],
                                of = of,
                                print_header = print_header,
                                other_info = other_info)

    return block_id

def collen(df):
    df = df.astype(str)
    q = '-->'

    df['cds_loc'] = df['start'].astype(str)+'..'+df['end'].astype(str)
    direction = len('dir')
    len_plen = df['plen'].astype(str).str.len().max()
    len_pid = df['pid'].str.len().max() if df['pid'].str.len().max() >= len('pid') else len('pid')
    len_type = df['type'].str.len().max() if df['type'].str.len().max() >= len('type') else len('type')
    len_gene = df['gene'].str.len().max() if df['gene'].str.len().max() >= len('gene') else len('gene')
    len_cds = df['cds_loc'].str.len().max()
    len_locus = df['locus'].str.len().max() if df['locus'].str.len().max() >= len('locus') else len('locus')

    len_modified = len('gi')
    len_product = 0

    collen_len = [int(x) for x in [len(q), len_cds, len('dir'), len_plen, int(len_pid), len_type,
                  len_gene, len_locus, len_modified, len_product]]
    return collen_len

if __name__ == '__main__':
    args = parse_cli()

    # Setting arguments
    Entrez.email = NcbiConfig['email'] if exists email in NcbiConfig else args.email
    if not Entrez.email or Entrez.email == '':
        print(f'No e-mail! Set one at configuration file ~/.rotifer/etc/db/ncbi.yaml or use option --email', file=sys.stderr)
        exit(1)
    accs = args.accession              # List of accessions
    above = int(args.above)            # Get n neighbors upstream
    below = int(args.below)            # Get n neighbors downstream
    verbose = args.verbose             # Verbose
    outformat = args.outformat.lower() # Output format (Table or gi2operon)
    distance = int(args.distance) * (-1)
    search = args.search
    if verbose:
        verbose_msg('Starting acc2operon')

    # Load NCBI API key
    api_key = ''
    if args.api_key:
        api_key = args.api_key
    elif "NCBI_API_KEY" in os.environ:
        api_key = os.environ["NCBI_API_KEY"]
    else:
        try:
            api_key = rcf.loadAPI(args.api_key)
            #verbose_msg(f'api: {api_key}')
            #exit()
        except:
            verbose_msg(f'No NCBI_API_KEY found. Please set your enviroment or use --api_key')
    if args.verbose: verbose_msg(f'NCBI_API_KEY: {api_key}')

    # Should remove this
    # Load previous accs search table
    df_accs = pd.DataFrame(columns = ['ids', 'acc', 'nucleotide', 'assembly'])

    header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene product organism classification'.split()
    # Include columns

    # Add other info to gi2operon output
    other_info = args.addinfo # Gi2operon header

    asms = args.assembly
    # Debbuging
    # Download asm and gbk to cache folder
    block_id = 0

    # Filter by assembly
    q_asm = ''
    if asms:
        # print(accs)
        asm = rcf.openread(asms)
        q_asm = f"AND assembly in ({', '.join(asm)})"

    # Filter by nucleotide
    q_nuc = ''
    # if nucs:
    #     nuc = rcf.openread(nucs)
    #     q_asm = f"AND nucleotide in ({', '.join(nuc)})"

    if search == 'local':
        try:
            click = clickhouse(table_name = 'genomes')
            conn = click.conn

            # This is a tuple nucleotide, nuc_asm
            s = time.time()
            acc_formated = ','.join(["'" + x +"'" for x in accs])
            nucleotide_in_db = conn.execute(f"select distinct nucleotide, nuc_asm FROM genomes where pid in ({acc_formated})").fetchall()
            e = time.time()

            tot_group = len(nucleotide_in_db)
            left_group = 0
            if verbose:
                verbose_msg('Total assemblies: {0}'.format(str(tot_group)))

            if args.progress:
                for nucleotide, nuc_asm in tqdm(nucleotide_in_db):
                    block_id = block_intervals(accs = accs, above = above, below = below, block_id = block_id, of = outformat,
                                         conn = conn,
                                         nucleotide = nucleotide,
                                         nuc_asm = nuc_asm,
                                         other_info = other_info,
                                               distance = distance
                                         )
                    left_group +=1
                if verbose:
                    verbose_msg('Done')

            else:
                for nucleotide, nuc_asm in nucleotide_in_db:
                    block_id = block_intervals(accs = accs, above = above, below = below, block_id = block_id, of = outformat,
                                         conn = conn,
                                         nucleotide = nucleotide,
                                         nuc_asm = nuc_asm,
                                         other_info = other_info,
                                         distance = distance
                                         )
                    left_group +=1
                    if verbose:
                        verbose_msg(f'Remaining assemblies: {tot_group-left_group}')

        except:
            pass
        finally:
            try:
                if args.fun:
                    sys.stderr.write(args.fun.fun)
            except:
                pass

    else: # www + insert into local clickhouse
        # Create a tmp table name
        table_name = ''.join(random.choices(string.ascii_letters, k=24))

        try:
            if verbose:
                verbose_msg('Collecting ipgs and gbks')

            click = clickhouse(table_name = 'genomes')
            conn = click.conn

            # What does this do?
            df_all, df_ann =  acc2asm2gbk(accs, verbose, df_accs, threads = int(args.threads),
                              api_key = api_key, progress = args.progress)

            # Show statistics
            if verbose:
                verbose_msg('IPG reports collected! Summary statistics:')
                verbose_msg(f'df_ann size {df_ann.memory_usage().sum()/(1024**3)}')
                verbose_msg(f'df_ann_selected number rows {df_ann.shape[0]}')

            # Save IPG dataframes
            if args.ipg:
                verbose_msg('Saving IPG reports to '+args.ipg+'.{all,ann}.a2o.tsv')
                df_all.to_csv(args.ipg+'.all.a2o.tsv', sep = '\t', index = None)
                df_ann.to_csv(args.ipg+'.ann.a2o.tsv', sep = '\t', index = None)

            # Save missing dataframes
            if args.missing and missing:
                verbose_msg('List of protein accessions without IPG report saved to '+args.missing+'.missing.a2o.tsv')
                missing = set([x for x in accs if x not in set(df_all.original)])
                missing.sort()
                mfh = open(args.missing+'.missing.a2o.tsv', 'w')
                mfh.write("\n".join(missing))
                mfh.close()

            # What assemblies are NOT in our clickhouse server? Answer: ls_gbks
            ls_gbks = df_all['assembly_modified'].unique()
            ls_gbks_parsed = rcf._flatten([list(x) for x in conn.execute('''SELECT DISTINCT nuc_asm from genomes''').fetchall()])
            ls_gbks = list(set([x for x in ls_gbks if x not in ls_gbks_parsed]))
            n = len(ls_gbks)//int(args.threads) if len(ls_gbks)//int(args.threads) > 0 else 1
            sub_ls_gbks = [ls_gbks[x:x+n] for x in range(0,len(ls_gbks), n)]

            # Download missing assemblies
            jobs = []
            for x in range(len(sub_ls_gbks)):
                p = Process(target = _multi_ft2df, args = (sub_ls_gbks[x], verbose))
                jobs.append(p)
                p.start()
            try:
                for p in jobs:
                    p.join()
            except KeyboardInterrupt:
                for p in jobs:
                    p.terminate()
                sys.exit(2)

            # Get unique nucleotides
            if verbose:
                verbose_msg(f'Collected {len(ls_gbks_parsed)} gbks')

            block_id = 0

            # df_all.to_csv('df_all.test', sep = '\t', index = None)
            df_all = df_all[df_all['acc'].notnull()]

            accs_in_db = list(df_all['acc'].unique())

            if verbose:
                verbose_msg(f'Ran accs_in_db ok!')

            s = time.time()

            nucleotide_in_db = []
            n_partition = 2000
            # nucleotide_in_db = conn.execute(f"select distinct nucleotide, nuc_asm FROM genomes where pid in ({acc_formated})").fetchall()
            # This is a tuple nucleotide, nuc_asm
            if verbose:
                verbose_msg('Collecting nucleotides from pids using clickhouse')
            if args.hit == 'best':
                nucs_found = ["'" + x + "'" for x in df_all['assembly_modified'].unique()]
                for x in range(0,len(accs_in_db), n_partition):
                    # acc_formated = ','.join(["'" + x +"'" for x in accs])
                    acc_formated = ','.join(["'" + z +"'" for z in accs_in_db[x:x+n_partition]])
                    if verbose:
                        verbose_msg(f'{accs_in_db[x:x+n_partition][0:10]}...')
                        verbose_msg(f'{acc_formated[0:30]}...')
                        verbose_msg(f'last {accs_in_db[x]}')

                    try:
                        _ = conn.execute(f"""SELECT DISTINCT nucleotide, nuc_asm
                                                            FROM genomes
                                                            WHERE nuc_asm IN ({','.join(nucs_found)})
                                                                AND pid IN ({acc_formated})"""
                                                    )

                        nucleotide_in_db.extend([(nuc, nuc_asm) for nuc, nuc_asm in _])
                    except:
                        verbose_msg('Error')
                        # print('\n'.join(accs_in_db[x:x+n_partition]))

            elif args.hit == 'any':
                for x in range(0,len(accs_in_db), n_partition):
                    acc_formated = ','.join(["'" + z +"'" for z in accs_in_db[x:x+n_partition]])
                    _ = conn.execute(f"""SELECT DISTINCT nucleotide, nuc_asm
                                                        FROM genomes
                                                        WHERE pid IN ({acc_formated})"""
                                                )

                    nucleotide_in_db.extend([(nuc, nuc_asm) for nuc, nuc_asm in _])

            nucleotide_in_db = set(nucleotide_in_db)
            if verbose:
                verbose_msg(f'Ran nucleotide_in_db ok!')

            e = time.time()

            # Reduce df_ann size
            df_ann = df_ann[df_ann['nucleotide'].isin([x for x,y in nucleotide_in_db])]

            tot_group = len(nucleotide_in_db)
            left_group = 0

            if verbose:
                verbose_msg('Total assemblies: {0}'.format(str(tot_group)))

            if args.progress:
                for nucleotide, nuc_asm in tqdm(nucleotide_in_db):

                    left_group +=1
                    try:
                        original = df_ann[df_ann['nucleotide'] == nucleotide][['acc', 'representative']]

                        block_id = block_intervals(accs = accs_in_db, above = above, below = below, block_id = block_id, of = outformat,
                                             conn = conn,
                                             nucleotide = nucleotide,
                                             nuc_asm = nuc_asm,
                                             other_info = other_info,
                                             original_acc = original,
                                             distance = distance
                                             )

                    except:
                        pass

            else:

                for nucleotide, nuc_asm in (nucleotide_in_db):
                    left_group +=1
                    try:
                        original = df_ann[df_ann['nucleotide'] == nucleotide][['accession','representative']]

                        block_id = block_intervals(accs = accs_in_db, above = above, below = below, block_id = block_id, of = outformat,
                                             conn = conn,
                                             nucleotide = nucleotide,
                                             nuc_asm = nuc_asm,
                                             other_info = other_info,
                                             original_acc = original,
                                             distance = distance
                                             )


                        if verbose:
                            verbose_msg(f'Remaining assemblies: {tot_group-left_group}')
                    except:
                        pass
                        # if verbose:
                        #     verbose_msg(f'Error {nucleotide}, {nuc_asm}')

            client.execute(f'drop table rotifer.{table_name}')

        except:
            try:
                client.execute(f'drop table rotifer.{table_name}')
            except:
                pass

        # Clean Up
        finally:
            try:
                client.execute(f'drop table rotifer.{table_name}')
            except:
                pass

            try:
                if verbose:
                    verbose_msg('Cleaning cache')
                for gb in os.listdir(GlobalConfig['cache']):
                    os.remove(os.path.join(gbks_folder, gb))

            except:
                pass
            try:
                if args.fun:
                    sys.stderr.write(args.fun.fun)
            except:
                pass