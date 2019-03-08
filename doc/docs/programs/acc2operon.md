## Acc2operon

Retrieving neighborhood with acc2operon

Command:
`acc2operon [-options] <protein_accession>`

_____________________

### 1. Simple Usage

The simplest form to get a gene neighborhood is using a protein accession.

```bash
acc2operon AAF37898.1
```

Also, acc2operon accepts a file containing a list of protein accession. Each accession must be in a new line.

```bash
cat proteins.acc
>> ALA74808.1
>> ANB85424.1
>> ASO61569.1
>> AMS29978.1
>> ASB67180.1
>> NP_355129.1
>> AWK11510.1
>> AMQ86316.2
>> NP_417993.1
>> AKP44116.1

acc2operon proteins.acc
```

### 2. Options

Valid options includes:

<center>

|Option|Default value|
|:----:|------------ |
|[--above/-a](#2.1 above)|3|
|[--below/-b](#below)|3|
|[--threads/-t](#threads)|5|
|[--outformat/-of](#outformat)|table|
|[--search_method/-sm](#search-method)|www|
|[--hit_method/-hm](#hit-method)|best|
|[--progress/-p](#progress)||
|[--api_key](#api-key)|None|
|[--header/-y](#header)|None|
|[--configdump](#configdump)|None|
|[--configfile](#configfile)|None|

</center>

#### 2.01 Above

To collect more genes above the query, activate the flag `--above` or `-a` and insert how many genes above the query, the default value is `3`.

`acc2operon -a 3 <protein_accession>`

#### 2.02 below

To collect more genes below the query, activate the flag `--below` or `-b` and insert how many genes below the query, the default value is `3`.

`acc2operon -b 3 <protein_accession>`

#### 2.03 threads


#### 2.04 outformat


#### 2.05 search method


#### 2.06 hit method


#### 2.07 progress


#### 2.08 api key


#### 2.09 header


#### 2.10 configdump


#### 2.11 configfile


### 3. Dev

Developing this program includes:

- Develop the rotifer.genome.database module
- Develop the rotifer.genome.data module

#### 3.1 Block definition
We discussed about what is the minimun information that represents a genome data (block). This feature is particular important to model a better database (current using clickhouse, but keep in mind other DB must be included like MS SQL and postgres).

The minimun columns in a dataframe that defines an object belonging to rotifer.genome.data is:

- nucleotide
- type
- start
- end
- strand
- feature_id
- block_id

#### 3.2 Things to implement

##### Splicing and origin

A known issue regards when a feature is composed by more than one location (See Bio.SeqFeature, BioPython module).

When such feature occurs, like in spliced genes and/or features overlapping the origin this leads to an incorrect result.

To overcome this issue a possible solution is:
The signature of a compound feature is the presence of at least one pair of features, when sorted by position (start, end), includes:

- A feature whose end coordinate is equal to the maximum value of any end coordinate for the corresponding nucleotide/contig. (Feature_end = end_max)
- The next feature start coordinate is equal to the minimum start of all features of the corresponding nucleotide/contig. (Feature_start = start_min)

The presence of these properties implies that any feature, spliced or not, runs over the origin of a circular replicon.

##### Annotation reference

Add a table to store generic annotation and references to external data for each feature. Keep in mind that each feature if identified by the feature_order and the nucletide accession. That applies even to the compound/spliced features, that will span several rows in our main table, one row per location.


|feature_order|nucleotide_id|Attribute| Database | Value|
|-------------|-------------|---------|----------|------|
|1            |NZ121321.1   | db_xref |GeneID    | 316  |
|1            |NZ121321.1   | db_xref |taxon    |457400  |


##### Minimize the number of SQL queries

Try to minimize the number of SQL queries needed to identify the nucleotide accessions and genome_order intervals.

##### Add a new table to database

We might, in the future, use a table that related pairs of features with a qualifier, like:

|     |     |         |
|-----|-----|---------|
|gene1| cds1| coded_by|
|gene1| mrna1| transcript|

##### How to represent a CompoundFeature in the dataframe?

A compound feature, or a gene with multiple isoforms.

Internally the rotifer.genome.data must know how many isoforms per gene, where are the all possible valid combinations. One way to go is to represent several isoforms as CompoundFeatures. Each gene may have several CompoundFeatures.

For example, the user may want to know how many isoforms and the isoform protein sequence.

How to represent such abstract data? Multiple lines, each line representing one isoform? If a user select one gene it will represent several isoforms?

How we will represent it inside the rotifer.genome.data? It will be a new pandas datatype? How to implement this class (inheritance from BioSeq) inside a pandas dataframe, it will be slow?



