#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import argparse
from multiprocessing import Process, Pipe
from collections import defaultdict
import rotifer.core.cli as corecli
import rotifer.table.table as tb
import yaml
import argcomplete

__version__ = 0.3
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser(description = 'Get common columns in different msa')

    parser.add('-m','--msa',
                        dest = 'msa',
                        default = [],
                        action = 'append',
                        helper = 'Input MSA (usage: -m msa1 -m msa2)')

    parser.add('-of', '--outformat',
                dest = 'outformat',
                        default = 'table',
                        helper = 'Select output format (yaml/table)')

    # argcomplete.autocomplete(parser)
    args = parser.parse_args()

    return args


def aln(ls, send_end):
    result = []
    for col,col2 in ls:
        tmp = pd.DataFrame({'t':df_t[col].values,'s':df_msa[col2]})
        tmp['count'] = np.where(tmp.t == tmp.s, 1 , 0)

        if tmp['count'].sum()  == tmp.shape[0]:
            result.append([col,col2, 1])

    send_end.send(result)

if __name__ == '__main__':
    args = parse_cli()
    df_t = tb.fasta2df(open(args.msa[0]).read().splitlines(), split_seq = True)
    df_msa = tb.fasta2df(open(args.msa[1]).read().splitlines(), split_seq = True)
    combs = [(x,y) for x in df_t.columns for y in df_msa.columns if x != 'Header' or x != 'ID' or y != 'Header' or y != 'ID']
    t = 5
    length = len(combs) //t if len(combs) //t >0 else 1
    sub_combs = [combs[x:x+length] for x in range(0, len(combs), length)]
    jobs = []
    pipe_list = []

    for x in range(len(sub_combs)):
        recv_end, send_end = Pipe(False)
        p = Process(target = aln, args = (sub_combs[x], send_end))
        pipe_list.append(recv_end)
        jobs.append(p)
        p.start()

    try:
        for p in jobs:
            p.join()
    except:
        for p in jobs:
            p.terminate()
        sys.exit(2)
    result_list = [x.recv() for x in pipe_list]
    if args.outformat == 'table':
        dc = defaultdict(list)
        for ls in result_list:
            for k,v,_ in ls:
                dc[k].append(v)
        print('MSA-1\tMSA-2')
        for k,v in dc.items():
            if k not in ['Header', 'ID']:
                print(str(k) + '\t' + '; '.join([str(x) for x in v]))
    if args.outformat == 'yaml':
        dc = defaultdict(list)
        for ls in result_list:
            for k,v, _ in ls:
                if k not in ['Header', 'ID']:
                    dc['postion'].append(v)
        print(yaml.dump(dc, default_flow_style=False))
