#Rotifer 0 rlineage

## NAME
rlineage - Convert NCBI lineage column

## SYNOPSIS
**rlineage** -ic [*column*] -xc [*column*] -c [*column*] -t [*taxonomy*] [*FILE*]
\# Simplest usage:
**rlineage** rneighbors.tsv > rneighbors.tax.tsv

\# Advanced usage: select taxonomic classification using -t/--taxonomy
**rlineage** -t superkingdom -t phylum -t class -t order rneighbors.tsv > rneighbors.tax.tsv

## DESCRIPTION
This program extract taxonomic unit from a NCBI lineage column

## AUTHORS
Gilberto Kaihami

## OPTIONS
    This section describes program options, aliases and default values
    **Program options**
    **table**
    Input table file. Pipe compatible
    default: <empty>

    **--include_column**, **-ic**
    


## Program options summary

|Option|alias|value|
|-----:|:----|-----|
|a|haha|str|


\# Simplest case: use default options
rlineage rneighbors.tsv > rneighbors.tax.txt
HAHAHAHA

Get taxonomic unit from NCBI lineage format

rlineage rneighbors.tsv > rneighbors.tax.txt

