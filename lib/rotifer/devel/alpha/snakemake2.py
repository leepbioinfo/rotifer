from rotifer.db import ncbi
import pandas as pd
input_file = snakemake.input["acc"]
output_file = snakemake.output["ipg_pickle"]

with open("accessions.txt") as f:
    pid = [line.strip() for line in f if line.strip()]

ic = ncbi.IPGCursor()
ipgtable = ic.fetchall(pid)
i1 = ipgtable.query('pid in @pid')
missing = list(set(pid) - set(i1.pid))
i2 = ipgtable.query('representative in @missing')
i2.pid = i2.representative
ipg = pd.concat([i1,i2]).drop_duplicates('pid')
ipg.to_pickle(output_file)
