runhmmscan() {
	target=$1
	name=$2
	input=$3
	database=$4
	
	if [ -f ${target}.phobius.tsv ]; then
		phobius=${target}.phobius.tsv
	else
		phobius=""
	fi

	cut -f1 -d " " ${input} \
	| parallel --pipe -N1 -j36 --recstart '>' hmmscan --cpu 1 ${database} - \
	> ${target}.${name}.hmmscan.out 2> ${target}.${name}.hmmscan.err
	
	hmmer2table ${target}.${name}.hmmscan.out > ${target}.${name}.hmmscan.tsv
	domain2architecture -e 0.0101 ${phobius} ${target}.${name}.hmmscan.tsv > ${target}.${name}.hmmscan.arch
	architecture2table ${target}.${name}.hmmscan.arch > ${target}.${name}.hmmscan.arch.tsv
	ln -sf ${target}.${name}.hmmscan.arch ${target}.${name}.scan.arch
	ln -sf ${target}.${name}.hmmscan.arch.tsv ${target}.${name}.scan.arch.tsv
}

