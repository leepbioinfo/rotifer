function dom_count() {
    arch=$1
    shift
    architecture2table $arch | tgroup -f '$F[2]=$F[1];1' -g 0 -g 1 -a 2=count -i 0..2 "$@"
}
