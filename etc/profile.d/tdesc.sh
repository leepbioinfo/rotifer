function tdesc() {
    sep=$2;
    if [ "$sep" == "" ]; then
        sep='\t';
    fi;
    head -n1 $1 | ttranspose -s $sep | tfilter -f 'unshift(@F,$H{c}++);1';
}
