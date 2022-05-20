#!/usr/bin/env python3

def hhmdb(
    df,
    esl_index_file,
    grouper="c80e3",
    redundancy_cluster="c80i70",
    cpu=10,
    prefix="hhmdb",
    size=3,
):
    '''
    Builds a3m and hhm database for the chosen cluster in the current directory
    '''

    import tempfile
    from subprocess import Popen, PIPE, STDOUT
    import os

    curr_dir = os.getcwd()
    if not os.path.exists(curr_dir+"/hhmdb/aln/"):
        os.makedirs(curr_dir+"/hhmdb/aln/")

    with tempfile.TemporaryDirectory() as tdf:
        for group_cluster in df[grouper].unique():
            if df[df[grouper] == group_cluster][redundancy_cluster].drop_duplicates().dropna().nunique() >= size :
                df[df[grouper] == group_cluster][
                redundancy_cluster
                ].drop_duplicates().dropna().to_csv(
                f"{tdf}/accs", index=None, header=None
                )
                Popen(
                f"esl-sfetch -f {esl_index_file} {tdf}/accs | mafft --anysymbol --maxiterate 1000 --localpair --thread {cpu} - > {tdf}/accs.fa",
                stdout=PIPE,
                shell=True,
                ).communicate()
                sequence(f"{tdf}/accs.fa").to_file(f"./hhmdb/aln/{group_cluster}.fa")
    Popen(f'for x in ./hhmdb/aln/*.fa; do y=$(echo $x|cut -f 4 -d "/" | sed "s/\.fa//"); echo \#$y|cat - $x > ./hhmdb/aln/$y.aln; rm $x;done', stdout=PIPE, shell=True).communicate()
    with tempfile.TemporaryDirectory() as td:
        Popen(f"ffindex_build -s {td}/{prefix}.msa.ffdata {td}/{prefix}.msa.ffindex ./hhmdb/aln/", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_apply {td}/{prefix}.msa.ffdata {td}/{prefix}.msa.ffindex -i {td}/{prefix}.a3m.ffindex -d {td}/{prefix}.a3m.ffdata -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_apply {td}/{prefix}.a3m.ffdata {td}/{prefix}.a3m.ffindex -i {td}/{prefix}.hhm.ffindex -d {td}/{prefix}.hhm.ffdata -- hhmake -i stdin -o stdout -v 0", stdout=PIPE, shell=True).communicate()
        Popen(f"cstranslate -f -x 0.3 -c 4 -I a3m -i {td}/{prefix}.a3m -o {td}/{prefix}.cs219", stdout=PIPE, shell=True).communicate()
        Popen(f" sort -k3 -n -r {td}/{prefix}.cs219.ffindex | cut -f1 > {td}/{prefix}.sorting.dat", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_order {td}/{prefix}.sorting.dat {td}/{prefix}.hhm.ffdata {td}/{prefix}.hhm.ffindex {td}/{prefix}.hhm.ordered.ffdata {td}/{prefix}.hhm.ordered.ffindex ", stdout=PIPE, shell=True).communicate()
        Popen(f"ffindex_order {td}/{prefix}.sorting.dat {td}/{prefix}.a3m.ffdata {td}/{prefix}.a3m.ffindex {td}/{prefix}.a3m.ordered.ffdata {td}/{prefix}.a3m.ordered.ffindex ", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.cs219.ffdata ./hhmdb/{prefix}_cs219.ffdata", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.a3m.ordered.ffdata ./hhmdb/{prefix}_a3m.ffdata", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.hhm.ordered.ffdata ./hhmdb/{prefix}_hhm.ffdata", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.cs219.ffindex ./hhmdb/{prefix}_cs219.ffindex", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.a3m.ordered.ffindex ./hhmdb/{prefix}_a3m.ffindex", stdout=PIPE, shell=True).communicate()
        Popen(f"mv {td}/{prefix}.hhm.ordered.ffindex ./hhmdb/{prefix}_hhm.ffindex", stdout=PIPE, shell=True).communicate()
        return()
