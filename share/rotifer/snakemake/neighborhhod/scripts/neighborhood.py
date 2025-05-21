from rotifer.db import ncbi
from rotifer.devel.alpha import gian_func as gf
import pandas as pd
input_file = snakemake.input["acc"]
ipg_pickle= snakemake.input["ipg"]
output_file = snakemake.output["out_file"]
batch_size = snakemake.config.get("batch_size", 100)

with open(input_file) as f:
    pid = [line.strip() for line in f if line.strip()]

gnc = ncbi.GeneNeighborhoodCursor(batch_size=batch_size)
ipg = pd.read_pickle(ipg_pickle)
neighborhood_df = gnc.fetchall(pid,ipgs=ipg.query('superkingdom !="Eukaryota"'))
euk = ipg.query('superkingdom=="Eukaryota"').rename({'position':'block_id', 'stop':'end', 'description':'product'}, axis=1).eval('query =1').eval('type ="CDS"')
euk = euk.reindex(columns=neighborhood_df.columns)
neighborhood_df = pd.concat([neighborhood_df, euk], ignore_index=True)
neighborhood_df = gf.single_tax2ndf(neighborhood_df)
neighborhood_df.to_pickle(output_file)
