function dom_per_prot() {
    arch=$1;
    shift;
    architecture2table $arch | tgroup -g 0 -a 1=count -i 0..1 "$@";
}