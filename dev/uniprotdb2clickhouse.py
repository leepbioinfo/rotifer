#!/usr/bin/env python3

import rotifer.core.cli as corecli
import rotifer.core.functions as rcf
from rotifer.core.log import log
from rotifer.alchemy.connect import clickhouse
from clickhouse_driver import Client
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
import random, string
from clickhouse_sqlalchemy import Table, make_session, get_declarative_base, types, engines
from multiprocessing import Process
import os

import rotifer.core.cli as corecli
from tqdm import tqdm
__version__ = 0.001
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser(description = 'Add uniprot data to clickhouse')

    parser.add( long_arg = '-f',
                short_arg = '--folder',
                dest = 'fi',
                nargs = None,
                default = None,
                arg_type = None,
                helper = 'Uniprotdb file',
                action = "store"
                )
    # Add another options here

    args = parser.parse_args()

    return args

def insert_data(res, client):

    client.execute(f'INSERT INTO rotifer.uniprotdb values', res)

def put_data(files):
    client = Client('localhost')

    uri = 'clickhouse://default:@localhost/rotifer'
    engine = create_engine(uri)
    session = make_session(engine)
    metadata = MetaData(bind=engine)
    metadata.reflect(bind = engine)
    for to_open in files:
        fi = open(to_open)

        for l in fi:
            data = l.rstrip().split('\t')
            data_dc = {'uniprotid': data[0],
                       'id_type': data[1],
                       'id': data[2]}

            fetch = engine.execute(f"""select DISTINCT uniprotid from uniprotdb where uniprotid = '{data[0]}'""").fetchone()
            if not fetch: # not empty
                insert_data([{'uniprotid': data[0],
                             'id_type': 'UniProtKB-AC',
                             'id': data[0]}], client)

            else:
                pass
            insert_data([data_dc], client)

if __name__ == '__main__':
    args = parse_cli()

    folder = args.fi
    all_files = []
    for fi in os.listdir(folder):
        path = os.path.join(folder,fi)
        all_files.append(path)

    threads = 8
    size = len(all_files)//threads if len(all_files) // threads > 0 else 1
    splited_files = [all_files[x:x+size] for x in range(0, len(all_files), size)]


    for x in range(len(all_files)):
        p = Process(target = put_data, args(all_files[x], ))
        p.start()
        jobs.append(p)

    try:
        for p in jobs:
            p.join()
    except KeyboardInterrupt:
        for p in jobs:
            p.terminate()
        sys.exit(2)

