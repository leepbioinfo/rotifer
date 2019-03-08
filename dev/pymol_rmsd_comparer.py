import pymol
import sys
from multiprocessing import Process

__version__ = 0.1
__authors__ = 'Gilberto Kaihami; Aureliano Guedes'

def corr(ls, structure_dict, methods = ['align', 'super', 'cealign']):
    try:
        for a,b in ls:
            res = []
            if 'align' in methods:
                aln= pymol.cmd.align(structure_dict[a], structure_dict[b])
                res.append(str(aln[0]))
            if 'super' in methods:
                spr = pymol.cmd.super(structure_dict[a], structure_dict[b])
                res.append(str(spr[0]))
            if 'cealign' in methods:
                ce = pymol.cmd.cealign(structure_dict[a], structure_dict[b])
                res.append(str(ce['RMSD']))
            print('\t'.join([a,b] +res))
            sys.stdout.flush()
    except:
        pass


try:
    f = sys.argv[3]
    threads = int(sys.argv[4])
    opts = sys.argv[5:]
    _ = open(f, "r").read().splitlines()
    chains = {i.split('_')[0]: i.split('_')[1] for i in _}

    got = dict()
    structure_dict= {}

    for i in chains.keys():
        a = i.lower() + chains[i].upper()
        if a not in got.keys():
            got[a] = 1
            structure_dict[a] = pymol.cmd.fetch(a)

    sub_dict = [(x,y) for x in structure_dict.keys() for y in  structure_dict.keys() if y != x]
    n = len(sub_dict)//threads if len(sub_dict)//threads > 0 else 1
    sub_process = [sub_dict[x:x+n] for x in range(0,len(sub_dict), n)]
    print('\t'.join('PDB1 PDB2 {0}'.format(' '.join(opts)).split() ))


    jobs = []
    for x in range(len(sub_process)):
        p = Process(target = corr, args = (sub_process[x], structure_dict, opts))
        jobs.append(p)
        p.start()
    try:
        for p in jobs:
            p.join()
    except KeyboardInterrupt:
        for p in jobs:
            p.terminate()
        sys.exit(2)

except:
    sys.exit()
