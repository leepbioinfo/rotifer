# IPython log file

import os
import pandas as pd
from glob import glob
from rotifer.devel.beta import sequence as rdbs

# Loading our alignments
class SequenceCollection():
    def __init__(
            self,
            basedir = "/projects/salmonella/work/st",
            suffix = "sto",
            informat = "stockholm",
        ):
        alnpaths = glob(f'{basedir}/alndb/*.{suffix}')
        alndb = []
        alignments = dict()
        for x in alnpaths:
            basename = os.path.basename(x).replace(f".{suffix}","")
            pdb = glob(os.path.join(basedir,"colabfold") + "/" + basename + "/*_rank_001_*.pdb")
            if len(pdb):
                pdb = pdb[0]
            else:
                pdb=None
            alndb.append((basename, x, pdb))
            alignments[basename] = rdbs.sequence(x,informat)
        alndb = pd.DataFrame(alndb, columns=['name','path','pdb'])
        alndb.index = alndb.name.tolist()
        self.df = alndb
        self.alignment = alignments

    def __setitem__(self,key,value):
        self.df.loc[value.name] = [value.name,None,None]
        self.alignment[key] = value

    def __getitem__(self,key):
        return self.alignment[key]

    def keys(self):
        return self.df.index.tolist()

    # Saving to HTML
    def tohtml(self, model, filename, sample=50, consensus=90):
        aln = self[model].copy()
        pdb = self.df.loc[model,"pdb"]

        keep = []
        if pdb:
            aln = aln.add_pdb(aln.df.iloc[0,0], pdb_file=pdb)
            keep = aln.df.iloc[0:2,0].tolist()

        if len(aln) > sample:
            keep = keep + aln.df[2:].id.sample(sample-2).tolist()
            aln = aln.filter(keep=keep) 

        aln.to_html(consensus, filename, annotations=aln.df.iloc[0,0], remove_gaps=aln.df.iloc[1,0], fixed_index=True)

