#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import threading
from threading import Thread
import rotifer.core.cli as corecli
import argparse
import subprocess
import pandas as pd
from multiprocessing import Pool, Process, Pipe
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Entrez
import time
import shutil
import numpy as np
import warnings
from os.path import expanduser
from datetime import datetime as dt
import argcomplete
from clickhouse_driver import Client
from sqlalchemy import create_engine, Column, MetaData, literal
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, Integer, String
import random, string
from clickhouse_sqlalchemy import Table, make_session, get_declarative_base, types, engines
import rotifer.core.functions as rcf
import uuid
warnings.filterwarnings('ignore')

__version__ = 0.43
__authors__ = 'Gilberto Kaihami, Aureliano Guedes, Gabriel Hueck'

########
# New features
########

# In version 0.42:
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
                        helper = 'Insert a valid NCBI API key or a valid config file username',
                        default = '')

    parser.add(long_arg = '--search_method',
               short_arg = '-sm',
               dest = 'search',
               helper = 'Select a search method, local or www',
               default = 'wwww')

    parser.add(short_arg = '-c',
                        long_arg = '--cache',
               dest = 'cache',
                        helper = 'A folder name that all gbs will be saved',
                        default = 'cache')

    parser.add(short_arg = '-a',
                        long_arg ='--above',
               dest = 'above',
                        helper = 'Rows above: maximum number of neighbors upstream of target loci (Default: 3)',
                        default = 3)

    parser.add(short_arg = '-b',
                        long_arg = '--below',
               dest = 'below',
                        helper = 'Rows below: maximum number of neighbors downstream of target loci (Default: 3)',
                        default = 3)

    parser.add(long_arg = '--missing',
               dest = 'missing',
                        helper = 'Print missing accs',
                        action = 'store_true')

    parser.add(short_arg = '-v', long_arg = '--verbose',
               dest = 'verbose',
                       action = 'store_true')

    parser.add(short_arg = '-t',long_arg = '--threads',
               dest = 'threads',
                        helper = 'Number of threads (Default: 5)',
                        default = 5)

    parser.add(short_arg = '-of',
               long_arg = '--outformat',
               dest = 'outformat',
                        helper = 'Output format (table/gi2operon) [Default: table]',
                        default = 'table')

    # parser.add('--clean',
    #                     help = 'Clean cache folder',
    #                     action = 'store_true')



    argcomplete.autocomplete(parser)

    args = parser.parse_args()
    return args

def verbose_msg(message = ''):
    now = dt.now().strftime('[%D %H:%M:%S]')
    if isinstance(message, list):
        message = ' '.join(message)
    sys.stderr.write('## {0} {1}\n'.format(now, message))
    sys.stderr.flush()

def nuc2gbk(gbks,db, verbose = False, api_key = ''):
    for nucleotide in gbks:
        if verbose:
            verbose_msg('Downloading {0} using efetch'.format(nucleotide))
        try:
            to_save = Entrez.efetch(db = 'nuccore',
                                    rettype = 'gbwithparts',
                                    retmode = 'text',
                                    id = nucleotide,
                                    api_key = api_key)

            out_handle = open(db + '/' + nucleotide + '.gbff', "w")
            out_handle.write(to_save.read())
            out_handle.close()
            to_save.close()
        except:
            if verbose:
                if isinstance(nucleotide, str):
                    sys.stderr.write(nucleotide+'\n')

def asm2ftp(gbks, db, verbose = False):
    for ele in gbks:
        if verbose:
            verbose_msg('Downloading {0} using ftp'.format(str(ele)))
        tries_download = 0
        while tries_download < 10:
            try:
                gca, f1, f2, f3 = ele[0:3], ele [4:7], ele[7:10], ele[10:13]
                sub2 = subprocess.Popen('''curl -l ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/* 2> /dev/null| \
            grep "{4}" 2> /dev/null'''.format(gca, f1, f2, f3, ele), shell = True, stdout = subprocess.PIPE)
                genome = sub2.communicate()[0].decode('utf-8').replace('\n', '')
                sub1 = subprocess.Popen('''curl -l ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/* 2> /dev/null |\
            grep "{4}" | parallel curl -l ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/{5}/ 2> /dev/null'''.format(gca, f1, f2,f3, ele, '{}'), shell = True, stdout = subprocess.PIPE)
                features = sub1.communicate()[0].decode('utf-8').splitlines()
                for feature in features:
                    file2save = db+'/'+ele
                    if file2save+'.gbff' not in os.listdir(db):
                        if '_genomic.gbff' in feature and '_cds' not in feature and '_rna_' not in feature:
                            download = 'curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/{4}/{5} --output {6}.gbff.gz 2> /dev/null'.format(gca, f1, f2, f3, genome, feature, file2save)
                            subprocess.Popen([download], shell = True).wait()
                            subprocess.Popen(['gunzip -f {0}.gbff.gz'.format(file2save)], shell = True).wait()
                break

            except:
                tries_download +=1
                time.sleep(0.3)


# def asm2epost:
#     pass
#
# def asm2efetch:
#
#     pass
#
# def asm2esearch:
#     pass

def acc2asm2gbk(accs, db, verbose, df_accs, threads = 5, api_key = '', table_name = ''):
    # Split acc in a list of size s

    #### TODO remove useless stuf here
    s = 1000
    accs2find = [acc for acc in accs if acc not in df_accs['acc'].values] # Ok
    subs = [accs2find[x:x+s] for x in range(0, len(accs2find), s)]        # Ok

    # Setting lists
    # Ok here
    got = []
    got_a = got.append
    gbks_folder = os.path.join(db, 'gbff')

    try:
        try:
            os.makedirs(gbks_folder)
        except:
            pass

        try:
            os.makedirs(accs_folder)

        except:
            pass

        for sub in subs:
            nucleotides = []
            nucleotides_a = nucleotides.append
            ipg_ids = []
            ipg_ids_a = ipg_ids.append

            asms = []
            asms_a = asms.append

            pids = []
            pids_a = pids.append

            ipg_ids_all = []
            ipg_ids_all_a = ipg_ids_all.append

            # Find ipgs using efetch
            tries = 0
            ## Try using epost if fail except using efetch

            try:
                request = Entrez.epost(db = 'protein', id = ",".join(sub))
                result = Entrez.read(request)
                webEnv = result['WebEnv']
                querykey = result['QueryKey']

                handle = Entrez.efetch(db='protein',
                                       rettype = 'ipg',
                                       retmode = 'text',
                                       webenv = webEnv,
                                       query_key = querykey,
                                       api_key = api_key)

                seq_handle = [x for x in handle.read().split('\n')[1:] if x != '']

                for seq in seq_handle:
                    e = seq.split('\t')
                    pids_a(e[6])
                    nucleotides_a(e[2])
                    asms_a(e[10])

                    if e[0] not in ipg_ids:
                        ipg_ids_a((e[0], e[6])) # Fix this
                    ipg_ids_all_a(e[0])

                handle.close()
                time.sleep(0.3)
                for ele in sub:
                    got.append(ele)

            except:
                while tries < 5:
                    sub = [x for x in sub if x not in got]
                    try:
                        for sub_acc in sub:
                            request = Entrez.epost(db = 'protein', id = sub_acc)
                            result = Entrez.read(request)
                            webEnv = result['WebEnv']
                            querykey = result['QueryKey']

                            handle = Entrez.efetch(db='protein',
                                                   rettype = 'ipg',
                                                   retmode = 'text',
                                                   webenv = webEnv,
                                                   query_key = querykey,
                                                   api_key = api_key)

                            seq_handle = [x for x in handle.read().split('\n')[1:] if x != '']
                            for seq in seq_handle:
                                e = seq.split('\t')
                                pids_a(e[6])
                                nucleotides_a(e[2])
                                asms_a(e[10])
                                if e[0] not in ipg_ids:
                                    ipg_ids_a((e[0], e[6])) # Fix this
                                ipg_ids_all_a(e[0])
                            handle.close()
                            got.append(sub_acc)
                            time.sleep(0.3)
                        break

                    except:
                        tries +=1
                        if tries == 4 and len(sub) > 0:
                            sub = sub[1:]
                            tries = 0
                        time.sleep(0.5)
                    if len(sub) == 0:
                        break

            # map
            acc_and_ids = {k:v for k,v in ipg_ids if v in got}
            ###

            tmp_df = pd.DataFrame({'ids': ipg_ids_all,'acc': pids, 'nucleotide':nucleotides, 'assembly':asms,
                                   'assembly_modified': ''})
            tmp_df['original'] = tmp_df.ids.map(acc_and_ids)
            tmp_df = tmp_df[['ids','acc', 'nucleotide', 'assembly', 'assembly_modified', 'original']]

            if (isinstance(tmp_df['original'].values[0], np.float64)):
                sys.stderr.write('No results\n')
                sys.exit()
            for g in got:
                check = tmp_df[tmp_df['original'] == g]['assembly']

                if len(list(set(check.values))) == 1 and list(set(check.values))[0] == '':
                    tmp_df.loc[check.index, 'assembly_modified'] = 'No_asm'

            df_accs = pd.concat([df_accs, tmp_df])
            df_accs = df_accs[['ids','acc', 'nucleotide', 'assembly', 'assembly_modified', 'original']]

            # map
            ### Check if no assembly
            ###################
            ####### Ate aqui ok
            # Drop ipgs

            df_accs_all = tmp_df[['ids', 'acc', 'assembly', 'nucleotide']].copy() # copy for all accs
            df_accs_all['assembly'] = df_accs_all['assembly'].replace('',np.nan)
            df_accs_all['assembly'].fillna(df_accs_all['nucleotide'], inplace = True)
            df_accs_all = df_accs_all[['ids', 'acc', 'assembly', 'nucleotide']]


            df_accs_all['dc'] = df_accs_all.apply(lambda x: {'acc': x.acc,
                                                             'nucleotide': x.nucleotide,
                                                             'assembly': x.assembly}, 1)
            # Insert acc
            client.execute(f'Insert into rotifer.{table_name} values', list(df_accs_all['dc'].values))


    #######

            first_ocurrence = df_accs[(df_accs['acc'].isin(sub)) & (df_accs['assembly_modified'] != 'No_asm')] # select
            first_ocurrence = first_ocurrence.drop_duplicates(subset = ['acc'], keep = 'first')

            asm2get = first_ocurrence['assembly'].values
    #        print(first_ocurrence)
    #        sys.exit()
            view_ids = first_ocurrence['ids'].values

            ##
            no_ocurrence = df_accs[~(df_accs['ids'].isin(view_ids))]
            no_ocurrence = no_ocurrence[no_ocurrence['assembly_modified'] != 'No_asm'].drop_duplicates(subset= ['ids', 'original'], keep = 'first')

            drop_nucleotide = df_accs[df_accs['assembly_modified'] == 'No_asm'].drop_duplicates(subset = ['ids', 'original'],keep = 'first')
            drop_nucleotide['assembly'] = drop_nucleotide['nucleotide'].values
            df_accs = df_accs[(df_accs['assembly'].isin(asm2get))]
            df_accs = pd.concat([df_accs, drop_nucleotide, no_ocurrence])
            df_accs['assembly'].replace('', 'remove_me',inplace=True)
            df_accs = df_accs[df_accs.assembly.str.contains("remove_me") == False]
            ###

            uri = 'clickhouse://default:@localhost/rotifer'
            engine = create_engine(uri)
            session = make_session(engine)
            metadata = MetaData(bind=engine)
            metadata.reflect(bind = engine)

            # GBKs in the DB
            got_gbks = rcf._flatten([list(x) for x in engine.execute('''
                                SELECT DISTINCT nuc_asm from genomes
                                ''').fetchall()])

            # GBKs in the cache folder
            got_gbks_folder = [x.replace('.gbff', '').replace('.tsv','') for x in os.listdir(gbks_folder)]

            #######

            find_gbks = [x for x in list(set(df_accs['assembly'])) if x not in got_gbks_folder and x not in drop_nucleotide['assembly'].values]
            find_gbks = [x for x in find_gbks if x not in got_gbks and x not in drop_nucleotide['assembly'].values]
            #######

            find_nuc = [x for x in list(set(drop_nucleotide['nucleotide'])) if x not in got_gbks_folder]
            find_nuc = [x for x in find_nuc if x not in got_gbks]
            #######

            # Find gbk using assembly
            jobs = []
            n = len(find_gbks)//threads if len(find_gbks)//threads > 0 else 1
            sub_gbks = [find_gbks[x:x+n] for x in range(0, len(find_gbks), n)]
            for x in range(len(sub_gbks)):
                p = Process(target = asm2ftp, args = (sub_gbks[x], gbks_folder))
                jobs.append(p)
                p.start()

            try:
                for p in jobs:
                    p.join()

            except KeyboardInterrupt:
                for p in jobs:
                    p.terminate()
                sys.exit(2)

            if find_nuc:
                nuc2gbk(find_nuc, gbks_folder, api_key = api_key)

    except:
        if verbose:
            verbose_msg('Droping table in acc2asm2gbk')
        client.execute(f'drop table rotifer.{table_name}')

    try:
        engine.close()
    except:
        pass
    return df_accs

def _multiple_find_acc(sub_fi, gbk, db,  send_end):
    results = []
    for f in sub_fi:
        df = pd.read_csv(os.path.join(db,'acc/'+f), sep = '\t')
        df = df[df['assembly'] == gbk]['acc'].values
        results.extend(df)
    send_end.send(results)

def features2df(gbk, df, db = '', pair_acc = '', verbose = False):
    '''
    Parameters:
    gbk:  assembly or nucleotide
    df:   empty dataframe
    db:   Database folder
    pair_acc
    '''
#    sys.stderr.write(gbk+'\n')
    # , accs = sub_accs, pair_acc = old_new_acc

    # Create engine
    uri = 'clickhouse://default:@localhost/rotifer'
    engine = create_engine(uri)
    session = make_session(engine)
    metadata = MetaData(bind=engine)
    metadata.reflect(bind = engine)

    # Check if gbk is present in the database
    if not engine.execute(f"select DISTINCT nuc_asm FROM genomes where assembly = '{gbk}'").fetchone():
        if verbose:
            verbose_msg(f'## loading {gbk}')

        s = time.time()
        fts = SeqIO.index(os.path.join(db,'gbff')+'/'+gbk+'.gbff', "genbank")
        e = time.time()
        if verbose:
            verbose_msg(f'## total time to load: {e-s}')
        df = pd.DataFrame()
        not_types = ['gene', 'source', 'repeat_region']

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

                for (idx, ft) in enumerate(rec.features):
                    if 'plasmid' in ft.qualifiers:
                        seq_type = 'plasmid'
                    if ft.type not in not_types:
                        feature_type = ft.type
                        sub_fts = ft.qualifiers
                        pid = sub_fts['protein_id'][0] if 'protein_id' in sub_fts.keys() else ''
                        if pid == '' and 'CDS' in ft.type:
                            feature_type = 'PSE'

                        locus_tag = sub_fts['locus_tag'][0] if 'locus_tag' in sub_fts.keys() else '.'
                        location = ft.location
                        start = str(location.start+1).replace('>', '').replace('<', '')
                        end = str(location.end).replace('>', '').replace('<', '')

                        try:
                            strand = str(location.strand)
                        except:
                            strand = '.'

                        gene = sub_fts['gene'][0] if 'gene' in sub_fts.keys() else '.'
                        plen = len(ft.qualifiers['translation'][0]) if 'translation' in ft.qualifiers.keys() else ''
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
                                    'nuc_asm': gbk })
                # write res here
                client = Client('localhost')
                # Insert
                s = time.time()

                if len(res) > 0:
                    if verbose:
                        verbose_msg(f'## Inserting into clickhouse DB {gbk}')
                    client.execute('Insert into rotifer.genomes values', res)
                    e = time.time()
                    if verbose:
                        verbose_msg(f'## Total time to insert {e-s}')

def _multi_ft2df(ls_gbks = [], db = '', verbose = False):
    for gbk in set(ls_gbks):
        if verbose:
            verbose_msg('Ft2df of {0}'.format(gbk))

        try:
            df = pd.DataFrame()
            features2df(gbk, df, db = db, verbose = verbose)
        except:
            if verbose:
                _ = 'Error ft2df {0}'.format(gbk)
                verbose_msg(str(_))



def cleanUp(folder):
    shutil.rmtree(folder)

###################################

def intervals(df, accs = '', above = 3, below = 3, distance = 0, block_id = 0, of = 'table'):
    df_CDS = df[df['type'] == 'CDS']
    df_CDS = df_CDS.reset_index()
    accs_indexes = df_CDS[df_CDS['pid'].isin(accs)].index
    intervals = []
    get_more_up = ''
    get_more_down = ''
    for value in accs_indexes:
        find_up = value - above
        find_down = value + below
        if intervals:
            if intervals[-1][1] - find_up+1 >= distance:
                if above+below+1 >= df_CDS.index.max():
                    intervals[len(intervals)-1] = [0, df_CDS.index.max()]
                else:
                    intervals[len(intervals) -1 ] = [intervals[len(intervals)-1][0], find_down]
            else:
                if above + below +1 >= df_CDS.index.max():
                    intervals.append([0, df_CDS.index.max()])
                else:
                    intervals.append([find_up, find_down])
        else:
            if above +below +1 >= df_CDS.index.max():
                intervals.append([0, df_CDS.index.max()])
            else:
                intervals.append([find_up, find_down])

    for up_index, down_index in intervals:
        block_id +=1
        sub_df = pd.DataFrame()
        get_more_up = ''
        get_more_down = ''
        if up_index < df_CDS.index.min():
            if df_CDS['topology'].values[0] != 'linear':
                get_more_up = df_CDS.index.max() + up_index
            up_index = df_CDS.index.min()
        if down_index > df_CDS.index.max():
            if df_CDS['topology'].values[0] != 'linear':
                get_more_down = down_index - df_CDS.index.max()
            down_index = df_CDS.index.max()
        if get_more_up:
            sub_df = df[(df['internal_id'] >=df_CDS.iloc[get_more_up+1].internal_id) &
                    (df['internal_id'] <=df_CDS.iloc[df_CDS.index.max()].internal_id)]
        if not sub_df.empty:
            _ = df[(df['internal_id'] >=df_CDS.iloc[up_index].internal_id) &
                    (df['internal_id'] <= df_CDS.iloc[down_index].internal_id)]
            sub_df = pd.concat([sub_df, _])

        else:
            sub_df = df[(df['internal_id'] >=df_CDS.iloc[up_index].internal_id) &
                    (df['internal_id'] <= df_CDS.iloc[down_index].internal_id)]
        if get_more_down:
            if not sub_df.empty:
                _ = df[(df['internal_id'] >= df_CDS.iloc[0].internal_id) &
                      (df['internal_id'] <= df_CDS.iloc[get_more_down-1].internal_id)]
                sub_df = pd.concat([sub_df, _])
            else:
                sub_df = df[(df['internal_id'] >= df_CDS.iloc[0].internal_id) &
                       (df['internal_id'] <= df_CDS.iloc[get_more_down-1].internal_id)]
        sub_df = sub_df.fillna('.')
        if of == 'gi2operon' or of == 'gi2operons':
            sub_df['plen'].astype(str)
            if '' in sub_df['plen'].values:
                sub_df.loc[sub_df[sub_df['plen'] == ''].index, 'plen'] = '.'

            if '' in sub_df['pid']:
                sub_df.loc[sub_df[sub_df['pid'] == ''].index, 'pid'] = '.'

            col_len = collen(sub_df)
            _ = [x for x in df['pid'].values if x in accs]
            print('ORGANISM', df.organism.values[0], 'accession no is',
                  df.nucleotide.values[0], 'Protein is', _[0])
            header = [".","cds","dir","len","pid","type","gene","locus","gi",'product']
            header = [header[i].ljust(col_len[i]) for i in range(len(header))]

            print('  '.join(header))
        for i, row in sub_df.iterrows():
            query = 1 if row.pid in accs else 0
            if of == 'table':
                asm2print = row.assembly
                print('\t'.join(([str(x) for x in [
                row.nucleotide, row.start, row.end, row.strand, block_id,
                query, row.pid, row.type, row.plen, row.locus, row.seq_type, asm2print,
                    row.gene, row['product'], row.organism, row.taxonomy
                                                  ]
                                  ]
                                 )
                                )
                      )
            if of == 'gi2operon' or of == 'gi2operons':
                query = '-->' if row.pid in accs else '.'
                direction = '+' if row.strand == 1 else '-'

                cds = str(row.start) +'..'+str(row.end)
                toprint = [str(x) for x in [query, cds, direction,
                                  row.plen, row.pid, row.type,
                                  row.gene, row.locus, ".", row['product']]]
                toprint = [toprint[i].ljust(col_len[i]) for i in range(len(toprint))]
                print('  '.join(toprint))
        if of == 'gi2operon' or of == 'gi2operons':
            print("---------------------------------------")

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

def makeDB(verbose = False):
    home = expanduser("~") #Get home folder platform independent
    db = os.path.join(home, '.cache/rotifer/acc2operon')
    try:
        os.makedirs(db)
        return db
    except:
        if verbose:
            sys.stderr.write('## Cache folder already exists at: {0}\n'.format(db))
        return db


if __name__ == '__main__':
    Entrez.email = "kaihami@usp.br"

    # Setting arguments
    args = parse_cli()
    accs = args.accession              # List of accessions
    above = int(args.above)            # Get n neighbors upstream
    below = int(args.below)            # Get n neighbors downstream
    verbose = args.verbose             # Verbose
    outformat = args.outformat.lower() # Output format (Table or gi2operon)
    db = makeDB(verbose)

    search = args.search
    if verbose:
        verbose_msg('Starting acc2operon')

    # Check API
    api_key = ''
    try:
        api_key = cf.loadAPI(args.api_key)

    except: pass

    if api_key == '':
        api_key = args.api_key

    if api_key == '':
        api_key = '935f49538aaa065207b4582a5cc4fcd59408'
        # sys.stderr.write('Insert valid NCBI API Key\n')

    # Should remove this
    # Load previous accs search table
    df_accs = pd.DataFrame(columns = ['ids', 'acc', 'nucleotide', 'assembly'])


    # Debbuging
    # Download asm and gbk to db folder
    block_id = 0
    if search == 'local':
        try:
            uri = 'clickhouse://default:@localhost/rotifer'
            engine = create_engine(uri)
            session = make_session(engine)
            metadata = MetaData(bind=engine)
            metadata.reflect(bind = engine)

            acc_formated = ','.join(["'" + x +"'" for x in accs])

            # Print header
            header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene product organism classification'.split()

            if outformat == 'table':
                print('\t'.join(header))
            # Get list of nucleotides

            nucleotide_in_db = rcf._flatten([list(x) for x in engine.execute(f"select DISTINCT nucleotide FROM genomes where pid in ({acc_formated})")])
            for nucleotide in nucleotide_in_db:
                df = pd.DataFrame(engine.execute(f"""select
                                                  internal_id,nucleotide,
                                                           assembly, topology, start, end,
                                                           strand,pid, type, plen, locus,
                                                           seq_type, gene, product,organism,
                                                           taxonomy
                                                 from genomes where nucleotide = '{nucleotide}' ORDER BY internal_id""").fetchall(),
                                                columns = ['internal_id','nucleotide',
                                                           'assembly', 'topology', 'start', 'end',
                                                           'strand','pid', 'type', 'plen', 'locus',
                                                           'seq_type', 'gene', 'product','organism',
                                                           'taxonomy'])
                df = df.drop_duplicates()

                block_id = intervals(df, accs, above, below, block_id = block_id, of = outformat)
        except:
            pass
        finally:
            try:
                if args.fun:
                    sys.stderr.write(args.fun.fun)
            except:
                pass

    else:
        # Create a tmp table name
        table_name = ''.join(random.choices(string.ascii_letters, k=24))

        try:

            if verbose:
                verbose_msg('Collecting ipgs and gbks')

            client = Client('localhost')

            client.execute(f'''CREATE TABLE rotifer.{table_name} (acc String, assembly String, nucleotide String)
                           ENGINE = MergeTree ORDER BY (acc, nucleotide)''')


            # Need to remove from db
            df_all = acc2asm2gbk(accs, db, verbose, df_accs, threads = int(args.threads),
                                 api_key = api_key, table_name = table_name)


            if verbose:
                verbose_msg('Collected ipgs and gbks')

            # Parse each gbk and print output
            #
            sub_df_all = df_all[df_all['original'].isin(accs)]

            # Parse gbks here

            ls_gbks = sub_df_all['assembly'].unique()
            uri = 'clickhouse://default:@localhost/rotifer'
            engine = create_engine(uri)
            session = make_session(engine)
            metadata = MetaData(bind=engine)
            metadata.reflect(bind = engine)
            ls_gbks_parsed = rcf._flatten([list(x) for x in engine.execute('''
                                SELECT DISTINCT nuc_asm from genomes
                                ''').fetchall()])

            jobs = []

            ls_gbks = list(set([x for x in ls_gbks if x not in ls_gbks_parsed]))
            n = len(ls_gbks)//int(args.threads) if len(ls_gbks)//int(args.threads) > 0 else 1
            sub_ls_gbks = [ls_gbks[x:x+n] for x in range(0,len(ls_gbks), n)]
            for x in range(len(sub_ls_gbks)):
                p = Process(target = _multi_ft2df, args = (sub_ls_gbks[x], db, verbose))
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

            sub_groups = sub_df_all['nucleotide'].unique()
            block_id = 0

            if args.missing:
                sys.stderr.write('## Missing acc\n')
                sys.stderr.write('\n'.join([x for x in accs if x not in df_all['acc'].values])+ '\n\n')
                sys.exit(0)
            header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene product organism classification'.split()

            if outformat == 'table':
                print('\t'.join(header))

            tot_group = len(sub_groups)
            left_group = 0
            if verbose:
                verbose_msg('Total assemblies: {0}'.format(str(tot_group)))

            for nucleotide in sub_groups:
                engine = create_engine(uri)
                session = make_session(engine)
                metadata = MetaData(bind=engine)
                metadata.reflect(bind = engine)
                #nucleotide
                df = pd.DataFrame(engine.execute(f"""select
                                                  internal_id,nucleotide,
                                                           assembly, topology, start, end,
                                                           strand,pid, type, plen, locus,
                                                           seq_type, gene, product,organism,
                                                           taxonomy
                                                 from genomes where nucleotide = '{nucleotide}' ORDER BY internal_id""").fetchall(),
                                                columns = ['internal_id','nucleotide',
                                                           'assembly', 'topology', 'start', 'end',
                                                           'strand','pid', 'type', 'plen', 'locus',
                                                           'seq_type', 'gene', 'product','organism',
                                                           'taxonomy'])

                df = df.drop_duplicates()
                find_accs = engine.execute(f"select acc from {table_name} where nucleotide = '{nucleotide}'").fetchall()

                accs_temp = rcf._flatten([list(x) for x in find_accs])
                block_id = intervals(df, accs_temp, above, below, block_id = block_id, of = outformat)

                left_group +=1

                if verbose:
                    verbose_msg(f'Remaining assemblies: {tot_group-left_group}')

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
                gbks_folder = os.path.join(db, 'gbff')
                for gb in os.listdir(gbks_folder):
                    os.remove(os.path.join(gbks_folder, gb))
            except:
                pass
            try:
                if args.fun:
                    sys.stderr.write(args.fun.fun)
            except:
                pass
