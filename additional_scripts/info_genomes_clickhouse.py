#!/usr/bin/env python3

from rotifer.alchemy.connect import clickhouse
import pandas as pd
import sys
from tabulate import tabulate

uri = 'clickhouse://default:@localhost/rotifer'
table_name = 'genomes'
db =  clickhouse(uri = uri,
                        table_name = table_name)
conn = db.conn
db_columns = db.table_columns()

# Select date
#
p = conn.execute('''select date, count(DISTINCT nuc_asm) from genomes GROUP BY date
                ORDER BY date''').fetchall()
df = pd.DataFrame(p, columns = ['date', 'genome_count'])


print('# Table')
size = conn.execute(f'''SELECT table, formatReadableSize(size) as size, rows FROM
        (SELECT table, sum(bytes) AS size, sum(rows) AS rows,
            min(min_date) AS min_date, max(max_date) AS max_date,
            (max_date - min_date) AS days
        FROM system.parts
        WHERE active     GROUP BY table     ORDER BY rows DESC )
    WHERE table = '{table_name}' ''').fetchall()[0]
print('Table name:         ', size[0])
print('Table size:         ', size[1])
print(f'Number of rows:     {size[2]:.2e}')
print(f'Distinct assmblies: {df["genome_count"].sum():.3e}')
print()

print('# Number of insertions per day')
print(tabulate(df, headers='keys', tablefmt='psql'))

