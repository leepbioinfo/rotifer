#Rotifer 0 rnexplorer

## NAME
rnexplorer - Annotate, clean, reformat rneighbors output

## SYNOPSIS
**rnexplorer** [*OPTIONS*] [*FILE*]

\# Simplest usage:

**rnexplorer** rneighbors.tsv -of gi2operon > rneighbors.gi2operon

## DESCRIPTION

This program helps you explore snapshots of gene neighborhooods created by  
rneighbors/acc2operon (Rotifer) and/or gi2operons (TASS)

## AUTHORS

Gilberto Kaihami; Robson Souza

## OPTIONS

This section describes program options, aliases and default values

**Program options**

**file**

Input table/gi2operon file. Pipe compatible

default: empty


**--header**, **-y**

Add more information to gi2operon header. The additional information must be a valid collumn  
in the table format

example: rnexplorer -y assembly -y classification

default: empty

**--outformat**, **-of**

Set output format. Acceptable output formats are:

    - table:     all neighborhood data formatted as a text table
    - gi2operon: all neighborhood data formatted as gi2operon
    - compact:   all neighborhood data formatted as compact

default: same as input

**--exclude\_by\_type**, **-xt**

Exclude individual features type in the format KEY:VALUE.  
The value could be an unique value (e.g. -xt type:CDS) or a file containing one column  
with values in each line.

example: -xt type:CDS
         -xt seq\_type:[FILE]
         -xt seq\_type:plasmid

**--include\_by\_type**, **-it**
Include individual features type in the format KEY:VALUE.  
The value could be an unique value (e.g. -it type:CDS) or a file containing one column  
with values in each line.

example: -it type:CDS
         -it seq\_type:[FILE]
         -it seq\_type:plasmid

**--annotation**, **-a**

Add annotation column to the output.  

Usage formats: [FILE] or [COL\_KEY], [COL\_KEY:FILE], or [COL\_KEY:FILE:COL\_NAME]  
The annotation file must contain at least two columns. The first column is the key and  
the second colum contains the annotation.  
There are three ways to use this flag:  
(1) Input just the annotation file it will map the first column against the pid column  
in the neighborhood table.  
(2) Include a column in the output table. The column can be selected by number, column  
name or interval (see exclude\_column)  
(3) You can pass the column to map and the file containing the annotation.  
(4) You can pass 3 parameters, the first is the reference column in the neighborhood  
table, the second is the annotation file, and the last parameter is the new column name.  

examples:

    (1) -a [FILE]    
    (2) -a pfam    
    (2) -a 0..3    
    (2) -a 1    
    (3) -a pid:[FILE]    
    (3) -a nucleotide:[FILE]    
    (4) -a pid:[FILE]:annotation\_col    

default: empty

**--padtable**, **-p**

Format gi2operon output

default: empty

**--exclude\_column**, **-xc**

Exclude a column in the output table  
Select a column by number (index starts from 0), column name or an interval  
Usage formats: col\_number, start..end, column name  

examples: -xc 0       (refers to the first column)
          -xc 0..4    (refers from the first to the fifth column)
          -xc product (refers to product column)

**--include\_column**, **-ic**

Include a column in the output table  
Select a columb by number (index starts from 0), column name or an interval  
Usage formats: col\_number, start..end, column name  

examples: -ic 0       (refers to the first column)
          -ic 0..4    (refers from the first to the fifth column)
          -ic product (refers to product column)

**--include\_column**, **-ic**

Include a column in the output table  
Select a columb by number (index starts from 0), column name or an interval  
Usage formats: col\_number, start..end, column name  

examples: -ic 0       (refers to the first column)
          -ic 0..4    (refers from the first to the fifth column)
          -ic product (refers to product column)

**--optionargs**, **-oa**

Advanced options.  
Select a architecture column for compact format in the formart 'arch':[COLUMN\_NAME] or 'arch'=[COLUMN\_NAME]  

example: rnexplorer -of compact -oa arch:pfam \<file\>  
         rnexplorer -of compact -oa arch=pfam -a pid:arch.tsv:pfam \<file\>  

## Program options summary

|Long name            |Aliases|Type   |
|--------------------|------|-----|
|--header             |-y     |list|
|--outformat          |-of    |string|
|--exclude\_by\_type  |-xt    |list|
|--include\_by\_type  |-it    |list|
|--annotation         |-a     |list |
|--padtable|-p|boolean|
|--exclude\_column    |-xc|list|
|--include\_column    |-it|list|
|--optionargs         |-oa|hash|
|--version|           |boolean|


