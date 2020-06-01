#!/bin/bash

### START OF TEMPLATE HEADER: do not change this section! 
function run_or_show() { # Do not use this function with pipes!
	if [ "$IS_TEST" == "1" ]; then
		echo "$@"
	else
		"$@"
	fi
}
#
# Rrsw arguments
IS_TEST=$1
SOURCE_DIR=$2
DEST_DIR=$3
export IWD=`pwd`
if [ "$IS_TEST" == "1" ]; then
	cd $DEST_DIR
fi
# LEEP environment
source /home/linuxbrew/.linuxbrew.sh
### END OF TEMPLATE HEADER

### START OF VARIABLE SECTION: you may change or add rows below

# Create taxid file for makeblastdb
#echo $(pwd) : $0 : $@

# Create taxonid maps
#if [ "$IS_TEST" == 1 ]; then
#	echo "\ls -1 /databases/taxonomy/accession2taxid/*prot.accession2taxid | xargs -n1 tail -n +2 | cut -f 2,3 | sed 's/\t/ /' > nr.taxid
#\ls -1 /databases/taxonomy/accession2taxid/*pdb.accession2taxid | xargs -n1 tail -n +2 | cut -f 2,3 | sed 's/\t/ /' > pdbaa.taxid
#\ls -1 /databases/taxonomy/accession2taxid/*.accession2taxid | grep -v prot | grep -v pdb | xargs -n1 tail -n +2 | cut -f 2,3 | sed 's/\t/ /' > nt.taxid
#"
#else
	#\ls -1 /databases/taxonomy/accession2taxid/*prot.accession2taxid | xargs -n1 tail -n +2 | cut -f 2,3 | sed 's/\t/ /' > nr.taxid
	#\ls -1 /databases/taxonomy/accession2taxid/*pdb.accession2taxid | xargs -n1 tail -n +2 | cut -f 2,3 | sed 's/\t/ /' > pdbaa.taxid
	#\ls -1 /databases/taxonomy/accession2taxid/*.accession2taxid | grep -v prot | grep -v pdb | xargs -n1 tail -n +2 | cut -f 2,3 | sed 's/\t/ /' > nt.taxid
#fi

# Prepare FASTA files for esl-sfetch
base=$(dirname $DEST_DIR)
if [ ! -d "$base/fadb"   ]; then run_or_show mkdir -p "$base/fadb"; fi
for f in $(\ls -1 db/FASTA/* | grep -v '\.md5$' | grep -v '\.ssi$' | grep -v '\.p..$')
do
	target=$(basename $f)
	taxid="${target}.taxid"
	dbtype="prot"
	if [ "$target" == "nt" ]; then DBTYPE="nucl"; fi
	run_or_show perl -i -pe 's/\-/X/go if (substr($_,0,1) ne ">")' $f
	run_or_show /home/linuxbrew/.linuxbrew/bin/esl-sfetch --index $f
	#run_or_show makeblastdb -hash_index -parse_seqids -blastdb_version 5 -dbtype prot -in $f -taxid_map $taxid -logfile makeblastdb.log
	#rm -f $taxid
        run_or_show ln -sf $DEST_DIR/db/FASTA/{$target,${target}.ssi} $base/fadb/
done

# Link databases!
if [ "$(\ls db/FASTA | grep -E '\.p..$')" != "" ]; then run_or_show ln -s db/FASTA/*.p?? .; fi
if [   -d "$base/freeze" ]; then run_or_show  ln -s $base/freeze .; fi

### END OF VARIABLE SECTION: your code ends here!

### START OF TEMPLATE FOOT: do not change this section! 
# Exit
cd $IWD
exit $?
### END OF TEMPLATE FOOT: any rows below will have no effect!

