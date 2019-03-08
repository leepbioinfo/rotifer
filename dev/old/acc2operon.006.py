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
from multiprocessing import Pool, Process
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
warnings.filterwarnings('ignore')

__version__ = 0.06
__authors__ = 'Gilberto Kaihami, Aureliano Guedes, Gabriel Hueck'

## TODO: get first valid nucleotide accession - protein acc - asm ==> Will allow to get PDB and uniprot!

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Get gene neighborhood',
                                     formatter_class = argparse.RawTextHelpFormatter) #, usage = msg())

    parser.add_argument('accession',
            help = "Input acession file",
        #    help=argparse.SUPPRESS,
    nargs = '*',
    action = corecli.action.autoload
    )
    parser.add_argument('-c',
                        '--cache',
                        help = 'A folder name that all gbs will be saved',
                        default = 'cache')
    parser.add_argument('-a',
                        '--above',
                        help = 'Rows above: maximum number of neighbors upstream of target loci (Default: 3)',
                        default = 3)
    parser.add_argument('-b',
                        '--below',
                        help = 'Rows below: maximum number of neighbors downstream of target loci (Default: 3)',
                        default = 3)
    parser.add_argument('--missing',
                        help = 'Print missing accs',
                        action = 'store_true')
    parser.add_argument('-v', '--verbose',
                       action = 'store_true')
    parser.add_argument('-t','--threads',
                        help = 'Number of threads (Default: 5)',
                        default = 5)
    parser.add_argument('-of','--outformat',
                        help = 'Output format (table/gi2operon) [Default: table]',
                        default = 'table')
    parser.add_argument('--clean',
                        help = 'Clean cache folder',
                        action = 'store_true')
    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'acc2operon',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = 'Get gene neighborhood'))
    args = parser.parse_args()
    return args


def nuc2gbk(gbks,db, verbose = False):
    for nucleotide in gbks:
        try:
            to_save = Entrez.efetch(db = 'nuccore', rettype = 'gbwithparts', retmode = 'text', id = nucleotide)
            out_handle = open(db + '/' + nucleotide + '.gbff', "w")
            out_handle.write(to_save.read())
            out_handle.close()
            to_save.close()
        except:
            if verbose:
                if isinstance(nucleotide, str):
                    sys.stderr.write(nucleotide+'\n')
def asm2ftp(gbks, db):
    for ele in gbks:
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
                    if '_genomic.gbff' in feature and '_cds' not in feature and '_rna_' not in feature:
                        file2save = db+'/'+ele
                        download = 'curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/{4}/{5} --output {6}.gbff.gz 2> /dev/null'.format(gca, f1, f2, f3, genome, feature, file2save)
                        subprocess.Popen([download], shell = True).wait()
                        subprocess.Popen(['gunzip {0}.gbff.gz'.format(file2save)], shell = True).wait()
                break
            except:
                tries_download +=1
                time.sleep(0.3)

def acc2asm2gbk(accs, db, verbose, df_accs, threads = 5):
    # Split acc in a list of size s
    s = 1000
    accs2find = [acc for acc in accs if acc not in df_accs['acc'].values] # Ok
    subs = [accs2find[x:x+s] for x in range(0, len(accs2find), s)]        # Ok
    if verbose:
        sys.stderr.write('## get asm\n')

    # Setting lists
    # Ok here
    got = []
    got_a = got.append
    gbks_folder = os.path.join(db,'gbff')
    try:
        os.makedirs(gbks_folder)
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

        tries = 0
        try:
            request = Entrez.epost(db = 'protein', id = ",".join(sub))
            result = Entrez.read(request)
            webEnv = result['WebEnv']
            querykey = result['QueryKey']
            handle = Entrez.efetch(db='protein', rettype = 'ipg', retmode = 'text', webenv = webEnv, query_key = querykey)
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
                        handle = Entrez.efetch(db='protein', rettype = 'ipg', retmode = 'text', webenv = webEnv, query_key = querykey)
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
#        df_accs.replace('N/A', np.nan, inplace = True)
#        df_accs.dropna(inplace = True)
        ###################
        ####### Ate aqui ok
        # Drop ipgs
        # CAJ84625.1 Problem here fix this
        first_ocurrence = df_accs[(df_accs['acc'].isin(sub)) & (df_accs['assembly_modified'] != 'No_asm')] # select
        asm2get = first_ocurrence['assembly'].values
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
        df_accs.to_csv(os.path.join(db, 'accs.tsv'), index = None, sep = '\t')
        ###

        got_gbks = [x.replace('.gbff', '') for x in os.listdir(gbks_folder)]
        find_gbks = [x for x in list(set(df_accs['assembly'])) if x not in got_gbks and x not in drop_nucleotide['assembly'].values]
        find_nuc = [x for x in list(set(drop_nucleotide['nucleotide']))]
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
            nuc2gbk(find_nuc, gbks_folder)
    got_gbks = [x.replace('.gbff', '') for x in os.listdir(gbks_folder)]
    df_accs = df_accs[df_accs['assembly'].isin(got_gbks)]
    df_accs.to_csv(os.path.join(db, 'accs.tsv'), index = None, sep = '\t')
    return df_accs

def features2df(gbk, df, db = '', accs = '', pair_acc = ''):
#    sys.stderr.write(gbk+'\n')
    fts = SeqIO.index(os.path.join(db,'gbff')+'/'+gbk+'.gbff', "genbank")
    types = ['CDS', 'tRNA']
    for k,rec in fts.items():
        idxs = []
        idxs_a = idxs.append
        pids = []
        pids_a = pids.append
        locus_tags = []
        locus_tags_a = locus_tags.append
        starts = []
        starts_a = starts.append
        ends = []
        ends_a = ends.append
        strands = []
        strands_a = strands.append
        genes = []
        genes_a = genes.append
        plens = []
        plens_a = plens.append
        products = []
        products_a = products.append
        feature_types = []
        feature_types_a = feature_types.append
        modified = []
        modified_a = modified.append
        seq_type = 'chromosome'
        gbk_id = k
        topology = rec.annotations['topology']
        organism = rec.annotations['organism']
        taxonomy = '; '.join(rec.annotations['taxonomy'])
        if rec.features:
            for ft in rec.features:
                if 'plasmid' in ft.qualifiers:
                    seq_type = 'plasmid'
            for ft in rec.features:
                if ft.type == 'CDS':
                    sub_fts = ft.qualifiers
                    pid = sub_fts['protein_id'][0] if 'protein_id' in sub_fts.keys() else ''
                    if pid in accs:

                        for (idx, ft) in enumerate(rec.features):
                            if ft.type in types:
                                feature_type = ft.type
                                sub_fts = ft.qualifiers
                                pid = sub_fts['protein_id'][0] if 'protein_id' in sub_fts.keys() else ''
                                if pid == '' and 'CDS' in ft.type:
                                    feature_type = 'PSE'
                                locus_tag = sub_fts['locus_tag'][0] if 'locus_tag' in sub_fts.keys() else '.'
                                location = ft.location
                                start = str(location.start+1)
                                end = str(location.end)
                                try:
                                    strand = str(location.strand)
                                except:
                                    strand = '.'
                                gene = sub_fts['gene'][0] if 'gene' in sub_fts.keys() else '.'
                                plen = len(ft.qualifiers['translation'][0]) if 'translation' in ft.qualifiers.keys() else ''
                                product = sub_fts['product'][0] if 'product' in sub_fts.keys() else ''
                                idxs_a(idx)
                                pids_a(pid)
                                locus_tags_a(locus_tag)
                                starts_a(start)
                                ends_a(end)
                                strands_a(strand)
                                genes_a(gene)
                                plens_a(plen)
                                products_a(product)
                                feature_types_a(feature_type)
                                space = len(pair_acc['original'].values[0])
                                original_acc = pair_acc[pair_acc['acc'] == pid]['original'].values[0] if pid in pair_acc['acc'].values and pid != pair_acc[pair_acc['acc'] == pid]['original'].values[0] else '.'
                                modified_a(original_acc)
                        _ = pd.DataFrame({'internal_id': idxs,
                                           'gbk': gbk_id,
                                           'topology': topology,
                                           'start':starts,
                                           'end': ends,
                                           'strand': strands,
                                           'pid': pids,
                                           'type': feature_types,
                                           'plen': plens,
                                           'locus': locus_tags,
                                           'seq_type': seq_type,
                                           'gene':genes,
                                           'product': products,
                                           'organism': organism,
                                           'taxonomy': taxonomy,
                                           'modified': modified
                        })
                        if df.empty:
                            df = _
                        else:
                            df = pd.concat([df, _])
                #            df = df.reset_index(drop=True)
                        break
    return (df.gbk.values, df)

def cleanUp(folder):
    shutil.rmtree(folder)

###################################

def intervals(df, above = 3, below = 3, distance = 0, block_id = 0, asm = '', of = 'table', pair_acc = ''):
    df_CDS = df[df['type'] == 'CDS']
    df_CDS = df_CDS.reset_index()
    accs = pair_acc['acc'].values
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
        if of == 'gi2operon' or of == 'gi2operons':
            sub_df['plen'].astype(str)
            if '' in sub_df['plen'].values:
                sub_df.loc[sub_df[sub_df['plen'] == ''].index, 'plen'] = '.'
            if '' in sub_df['pid']:
                sub_df.loc[sub_df[sub_df['pid'] == ''].index, 'pid'] = '.'
            col_len = collen(sub_df)
            _ = [x for x in df['pid'].values if x in accs]
            print('ORGANISM', df.organism.values[0], 'accession no is',
                  df.gbk.values[0], 'Protein is', _[0])
            header = [".","cds","dir","len","pid","type","gene","locus","gi",'product']
            header = [header[i].ljust(col_len[i]) for i in range(len(header))]
            print('  '.join(header))
        for i, row in sub_df.iterrows():
            query = 1 if row.pid in accs else 0
            original_acc = pair_acc[pair_acc['acc'] == row.pid] if row.pid in pair_acc['acc'].values else '.'
#            print(original_acc)
            if of == 'table':
                asm2print = asm if asm != row.gbk else '.'
                print('\t'.join(([str(x) for x in [
                row.gbk, row.start, row.end, row.strand, block_id,
                query, row.pid, row.type, row.plen, row.locus, row.seq_type, asm2print,
                    row.gene, row.modified.strip(), row['product'], row.organism, row.taxonomy
                                                  ]
                                  ]
                                 )
                                )
                      )
            if of == 'gi2operon' or of == 'gi2operons':
                query = '-->' if row.pid in accs else '.'
                direction = '+' if row.strand == '1' else '-'

                cds = str(row.start) +'..'+str(row.end)
                toprint = [str(x) for x in [query, cds, direction,
                                  row.plen, row.pid, row.type,
                                  row.gene, row.locus, row.modified, row['product']]]
                toprint = [toprint[i].ljust(col_len[i]) for i in range(len(toprint))]
                print('  '.join(toprint))
        if of == 'gi2operon' or of == 'gi2operons':
            print("---------------------------------------")

    return block_id
def collen(df):
    df['start'].astype(str)
    df['end'].astype(str)
    df['plen'].astype(str)
    q = '-->'
    df['cds_loc'] = df['start'].astype(str)+'..'+df['end'].astype(str)
    direction = len('dir')
    len_plen = df['plen'].astype(str).str.len().max()
    len_pid = df['pid'].str.len().max() if df['pid'].str.len().max() >= len('pid') else len('pid')
    len_type = df['type'].str.len().max() if df['type'].str.len().max() >= len('type') else len('type')
    len_gene = df['gene'].str.len().max() if df['gene'].str.len().max() >= len('gene') else len('gene')
    len_cds = df['cds_loc'].str.len().max()
    len_locus = df['locus'].str.len().max() if df['locus'].str.len().max() >= len('locus') else len('locus')

    len_modified = df['modified'].str.len().max() if df['modified'].str.len().max() >= len('gi') else len('gi')
    len_product = 0
    collen_len = [len(q), len_cds, len('dir'), len_plen, len_pid, len_type,
                  len_gene, len_locus, len_modified, len_product]
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

    # Load previous accs search table
    acc2load = os.path.join(db, 'accs.tsv')
    if os.path.isfile(acc2load):
        df_accs = pd.read_csv(acc2load, sep = '\t')
    else:
        df_accs = pd.DataFrame(columns = ['ids', 'acc', 'nucleotide', 'assembly'])

    # Debbuging
    # Download asm and gbk to db folder
    df_all = acc2asm2gbk(accs, db, verbose, df_accs, threads = int(args.threads))
    if verbose:
        sys.stderr.write('Finished acc2asm2gbk')
#    print(df_all)
    # Parse each gbk and print output
#    print('parsing gbk')
    sub_df_all = df_all[df_all['original'].isin(accs)]
#    print(sub_df_all)
    sub_groups = sub_df_all.groupby(by = ['assembly'])
    block_id = 0

    if args.missing:
        sys.stderr.write('## Missing acc\n')
        sys.stderr.write('\n'.join([x for x in accs if x not in df_all['acc'].values])+ '\n\n')
        sys.exit(0)
    header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene modified product organism classification'.split()
    if outformat == 'table':
        print('\t'.join(header))
    for i,val in sub_groups:
        gbk = val['assembly'].values[0]
        sub_accs = val['acc'].values
        old_new_acc = val
        df = pd.DataFrame()
        [nucleotide,df] = features2df(gbk, df, db = db, accs = sub_accs, pair_acc = old_new_acc)
        for nucleotide in set(nucleotide):
            block_id = intervals(df[df['gbk'] == nucleotide], above, below, asm = gbk, block_id = block_id, of = outformat, pair_acc = old_new_acc)
