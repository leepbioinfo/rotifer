# Table genome
create table genome ( internal_id Int64, nucleotide String, assembly String, topology String, start String, end String, strand String, pid String, type String, plen String, locus String, seq_type String, gene String, product String, organism String, taxonomy String, nuc_asm String) ENGINE = MergeTree ORDER BY (nuc_asm, pid, nucleotide, locus)
