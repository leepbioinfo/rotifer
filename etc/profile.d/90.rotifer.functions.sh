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

# List table columns
function tdesc() {
    sep=$2;
    if [ "$sep" == "" ]; then
        sep='\t';
    fi;
    head -n1 $1 | ttranspose -s $sep | tfilter -f 'unshift(@F,$H{c}++);1'
}

