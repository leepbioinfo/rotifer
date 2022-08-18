#!/bin/bash

set -e # Stop after first error

alndir=$1
prefix=$3

# Build initial database
ffindex_build -s $prefix.msa.ffdata $prefix.msa.ffindex $alndir/

# Add consensus
ffindex_apply $prefix.msa.ffdata $prefix.msa.ffindex \
 -i $prefix.a3m.ffindex -d $prefix.a3m.ffdata \
 -- hhconsensus -M 50 -maxres 65535 -i soutdirin -oa3m soutdirout -v 0

# Build HMM
ffindex_apply $prefix.a3m.ffdata $prefix.a3m.ffindex \
 -i $prefix.hhm.ffindex -d $prefix.hhm.ffdata \
 -- hhmake -i soutdirin -o soutdirout -v 0

# Calculate statistics
cstranslate -f -x 0.3 -c 4 -I a3m -i $prefix.a3m -o $prefix.cs219

# Reorder
sort -k3 -n -r $prefix.cs219.ffindex | cut -f1 > $prefix.sorting.dat
ffindex_order $prefix.sorting.dat $prefix.hhm.ffdata $prefix.hhm.ffindex \
 $prefix.hhm.ordered.ffdata $prefix.hhm.ordered.ffindex
ffindex_order $prefix.sorting.dat $prefix.a3m.ffdata $prefix.a3m.ffindex \
 $prefix.a3m.ordered.ffdata $prefix.a3m.ordered.ffindex

exit $?
