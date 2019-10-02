# Some Rotifer mantras turned into functions

# Get the number of domains for each protein (repetitions count as separate entries)
function dom_per_prot() {
    arch=$1
    shift
    architecture2table $arch | tgroup -g 0 -a 1=count -i 0..1 "$@"
}

# Get the number of UNIQUE domains for each protein (repetitions count only once)
function distinct_dom_per_prot() {
    arch=$1
    shift
    architecture2table $arch | tgroup -g 0 -a 1=count:distinct -i 0..1 "$@"
}

# Count domains
function dom_count() {
    arch=$1
    shift
    architecture2table $arch | tgroup -f '$F[2]=$F[1];1' -g 0 -g 1 -a 2=count -i 0..2 "$@"
}

function acc2pfam() {
    SEQ=$(efetch -db protein -format fasta -id $1)
    (
        echo "$SEQ" | hmmscan --cpu 4 /databases/pfam/Pfam - | hmmer2table -c;
        echo "$SEQ" | phobius 2> /dev/null | phobius2table -e 0.0101
    ) \
    | domain2architecture -e 0.0101 \
    | architecture2table \
    | padtable
}

function acc2profiledb() {
    (
        efetch -db protein -format fasta -id $1 | rpsblast -db ~/data/rpsdb/allprofiles | blast2table -c profiledb -s;
        efetch -db protein -format fasta -id $1 | phobius 2> /dev/null | phobius2table -e 0.0101
    ) \
    | domain2architecture -e 0.0101 \
    | architecture2table \
    | padtable
}

# List table columns
function tdesc() {
    sep=$2;
    if [ "$sep" == "" ]; then
        sep='\t';
    fi;
    head -n1 $1 | ttranspose -s $sep | tfilter -f 'unshift(@F,$H{c}++);1'
}

# Sort alignment by tree
function aln_order_by_tree() {
	aln=$1
	tree=$2
	tmp=$(mktemp -p /tmp aln_order_by_tree.XXXXXX)
	treeutil -of leaves $tree > $tmp
	aln2seqrows -r '\t' $aln | tjoin -r '\n' -i1 1 -i2 1 -1 1 -f1 '$F[1] ne "ID"' -f '$F[0]=">$F[0]";1' -p $tmp -
	rm -f $tmp
}

# Export all funcions
for f in $(grep ^function 90.rotifer.functions.sh | cut -f 2  -d " " | cut -f 1 -d '(')
do
	export -f $f
done
