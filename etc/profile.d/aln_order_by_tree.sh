function aln_order_by_tree() {
	aln=$1
	tree=$2
	tmp=$(mktemp -p /tmp aln_order_by_tree.XXXXXX)
	treeutil -of leaves $tree > $tmp
	aln2seqrows -r '\t' $aln | tjoin -r '\n' -i1 1 -i2 1 -1 1 -f1 '$F[1] ne "ID"' -f '$F[0]=">$F[0]";1' -p $tmp -
	rm -f $tmp
}
