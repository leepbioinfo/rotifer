function acc2pfam() {
    SEQ=$(efetch -db protein -format fasta -id $1);
    ( echo "$SEQ" | hmmscan --cpu 4 /databases/pfam/Pfam-A.hmm - | hmmer2table -c model=version; echo "$SEQ" | phobius 2> /dev/null | phobius2table -e 0.0101 ) | domain2architecture -e 0.0101 | architecture2table | padtable;
}
