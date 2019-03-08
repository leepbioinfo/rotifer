#!/usr/bin/env python3
from rotifer.alchemy.connect import clickhouse

click = clickhouse(table_name = 'genomes')
conn = click.conn


s = conn.execute('''
Select distinct(nuc_asm) from genomes2
''')

for e in s:
    if e[0] != '.':
        print(e[0])

