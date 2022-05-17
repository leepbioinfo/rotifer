#!/usr/bin/env python3

def cluster2hhdb(
    group_cluster,
    df,
    esl_index_file,
    grouper="c80e3",
    redundancy_cluster="c80i70",
    fast=False,
    cpu=8,
    view=False,
):
    '''
    Builds a3m and hhm database for the chosen cluster in the current directory
    '''

    import tempfile
    from subprocess import Popen, PIPE, STDOUT
    import os

    if not os.path.isdir("./hhdb/"):
        os.mkdir("hhdb")
        os.mkdir("./hhdb/"+group_cluster)
    with tempfile.TemporaryDirectory() as tmpdirname:
        df[df[grouper] == group_cluster][redundancy_cluster].drop_duplicates().dropna().to_csv(f"{tmpdirname}/accs", index=None, header=None)
        Popen(f"esl-sfetch -f {esl_index_file} {tmpdirname}/accs > {tmpdirname}/accs.fa", stdout=PIPE, shell=True).communicate()
        b = sequence(f"{tmpdirname}/accs.fa").realign(fast=fast, cpu=cpu)
    b.to_file(f"./hhdb/{group_cluster}/{group_cluster}.fa")
    Popen(f"ffindex_build ./hhdb/{group_cluster}/{group_cluster}.msa.ffdata ./hhdb/{group_cluster}/{group_cluster}.msa.ffindex ./hhdb/{group_cluster}/{group_cluster}.fa", stdout=PIPE, shell=True).communicate()
    Popen(f"ffindex_apply ./hhdb/{group_cluster}/{group_cluster}.msa.ffdata ./hhdb/{group_cluster}/{group_cluster}.msa.ffindex -i ./hhdb/{group_cluster}/{group_cluster}.a3m.ffindex -d ./hhdb/{group_cluster}/{group_cluster}.a3m.ffdata -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0", stdout=PIPE, shell=True).communicate()
    Popen(f"ffindex_apply ./hhdb/{group_cluster}/{group_cluster}.a3m.ffdata ./hhdb/{group_cluster}/{group_cluster}.a3m.ffindex -i ./hhdb/{group_cluster}/{group_cluster}.hhm.ffindex -d ./hhdb/{group_cluster}/{group_cluster}.hhm.ffdata -- hhmake -i stdin -o stdout -v 0", stdout=PIPE, shell=True).communicate()
    Popen(f"cstranslate -f -x 0.3 -c 4 -I a3m -i ./hhdb/{group_cluster}/{group_cluster}.a3m -o ./hhdb/{group_cluster}/{group_cluster}.cs219", stdout=PIPE, shell=True).communicate()
    Popen(f" sort -k3 -n -r ./hhdb/{group_cluster}/{group_cluster}.cs219.ffindex | cut -f1 > ./hhdb/{group_cluster}/{group_cluster}.sorting.dat", stdout=PIPE, shell=True).communicate()
    Popen(f"ffindex_order ./hhdb/{group_cluster}/{group_cluster}.sorting.dat ./hhdb/{group_cluster}/{group_cluster}.hhm.ffdata ./hhdb/{group_cluster}/{group_cluster}.hhm.ffindex ./hhdb/{group_cluster}/{group_cluster}.hhm.ordered.ffdata ./hhdb/{group_cluster}/{group_cluster}.hhm.ordered.ffindex ", stdout=PIPE, shell=True).communicate()
    Popen(f"ffindex_order ./hhdb/{group_cluster}/{group_cluster}.sorting.dat ./hhdb/{group_cluster}/{group_cluster}.a3m.ffdata ./hhdb/{group_cluster}/{group_cluster}.a3m.ffindex ./hhdb/{group_cluster}/{group_cluster}.a3m.ordered.ffdata ./hhdb/{group_cluster}/{group_cluster}.a3m.ordered.ffindex ", stdout=PIPE, shell=True).communicate()
    Popen(f"mv ./hhdb/{group_cluster}/{group_cluster}.a3m.ordered.ffdata ./hhdb/{group_cluster}/{group_cluster}.a3m.ffdata", stdout=PIPE, shell=True).communicate()
    Popen(f"mv ./hhdb/{group_cluster}/{group_cluster}.hhm.ordered.ffdata ./hhdb/{group_cluster}/{group_cluster}.hhm.ffdata", stdout=PIPE, shell=True).communicate()
    Popen(f"mv ./hhdb/{group_cluster}/{group_cluster}.a3m.ordered.ffindex ./hhdb/{group_cluster}/{group_cluster}.a3m.ffindex", stdout=PIPE, shell=True).communicate()
    Popen(f"mv ./hhdb/{group_cluster}/{group_cluster}.hhm.ordered.ffindex ./hhdb/{group_cluster}/{group_cluster}.hhm.ffindex", stdout=PIPE, shell=True).communicate()
    with open(f'./hhdb/{group_cluster}/{group_cluster}.hhm.ffdata') as f:
        hhm = f.read()
    if view:
        from IPython.core.page import page
        page(hhm)
    return(b, hhm)
