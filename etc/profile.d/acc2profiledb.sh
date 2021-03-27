function acc2profiledb() {
    (
        efetch -db protein -format fasta -id $1 | rpsblast -db $ROTIFER_DATA/rpsdb/allprofiles | blast2table -c profiledb -s;
        efetch -db protein -format fasta -id $1 | phobius 2> /dev/null | phobius2table -e 0.0101
    ) \
    | domain2architecture -e 0.0101 \
    | architecture2table \
    | padtable;
}
