#!/bin/bash
# File: /projects/salmonella/work/20210212/source1.sh

# Basename
target=$1
input=$2

# Detecting sequence signatures in the new proteins

# Phobius
cat ${input} | parallel -N1 -j10 --pipe --recstart ">" phobius > ${target}.phobius.out 2> ${target}.phobius.err
phobius2table -e 0.0101 ${target}.phobius.out > ${target}.phobius.tsv

# Pfam: hmmscan
cat ${target}.fa \
| parallel -N1 -j10 --pipe --recstart ">" hmmscan /databases/pfam/Pfam-A.hmm - \
> ${target}.pfam.hmmscan.out \
2> ${target}.pfam.hmmscan.err
hmmer2table ${target}.pfam.hmmscan.out > ${target}.pfam.hmmscan.tsv
domain2architecture -e 0.0101 ${target}.phobius.tsv ${target}.pfam.hmmscan.tsv > ${target}.pfam.hmmscan.arch
architecture2table ${target}.pfam.hmmscan.arch > ${target}.pfam.hmmscan.arch.tsv

# Aravind: hmmscan
cat ${target}.fa \
| parallel --pipe -N1 -j36 --recstart '>' hmmscan --cpu 1 /databases/profiledb/hmmer/aravinddb - \
> ${target}.aravind.hmmscan.out \
2> ${target}.aravind.hmmscan.err
hmmer2table ${target}.aravind.hmmscan.out > ${target}.aravind.hmmscan.tsv

# Aravind: rpsblast
cat ${target}.fa \
| parallel -N1 -j36 --pipe --recstart '>' rpsblast -db /databases/profiledb/rpsdb/aravinddb \
> ${target}.aravind.rpsblast.out \
2> ${target}.aravind.rpsblast.err
cat ${target}.aravind.rpsblast.out | parallel -N1 -j36 --pipe --recstart RPSBLAST blast2table -s > ${target}.aravind.rpsblast.tsv

# Aravind: domain2architecture
(cat ${target}.phobius.tsv; awk '$5<=0.1{print}' ${target}.aravind.rpsblast.tsv) \
| domain2architecture -e 0.0101 > ${target}.aravind.scan.arch
architecture2table ${target}.aravind.scan.arch > ${target}.aravind.scan.arch.tsv

