# IPython log file

import os
import pandas as pd
from glob import glob
from rotifer.devel.beta import sequence as rdbs
config = {
        "informat": {
            "*.sto": "stockholm",
            "*.fa": "fasta",
            "*.fasta": "fasta",
        }
}

# Loading our alignments
class SequenceCollection():
    def __init__(
            self,
            basedir = "/projects/salmonella/work/st/alndb",
            pattern = "*.sto",
            other_files = {
                "pdb": {
                    "basedir": "/projects/salmonella/work/st/colabfold",
                    "pattern": "/*_rank_001_*.pdb",
                }
            }
        ):
        if isinstance(pattern,str):
            pattern = [pattern]

        # Find alignments
        alndb = []
        alignments = dict()
        for extension in pattern:
            if "/**/" in extension:
                recursive = True
            else:
                recursive = False

            # Load data
            for x in glob(f'{basedir}/{extension}', recursive=recursive):
                # Cleanup target name
                if extension[0] == "*":
                    suffix = extension[1:]
                else:
                    suffix = extension
                basename = os.path.basename(x).replace(f"{suffix}","")
                alignments[basename] = rdbs.sequence(x,config["informat"][extension])

                # Load other paths
                otherfiles = []
                othernames = []
                for colname, other in other_files.items():
                    files = glob(os.path.join(other["basedir"], basename + other["pattern"]))
                    if len(files) == 0:
                        files=None
                    elif len(files) == 1:
                        files = files[0]
                    otherfiles.append(files)
                    othernames.append(colname)

                # Add data to stack
                alndb.append((basename, x, files))

        # Create internal dataframe for paths
        alndb = pd.DataFrame(alndb, columns=['name','path',*othernames])
        alndb.index = alndb.name.tolist()

        # Store data in object
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
    def tohtml(self, model, filename=None, sample=50, consensus=90):
        aln = self[model].copy()

        keep = []
        if "pdb" in self.df.columns:
            pdb = self.df.loc[model,"pdb"]
            if isinstance(pdb,list):
                pdb = pdb[0]
            aln = aln.add_pdb(aln.df.iloc[0,0], pdb_file=pdb)
            keep = aln.df.iloc[0:2,0].tolist()

        # Filter redundancy
        aln.add_cluster(coverage=0, identity=70, inplace=True)
        aln = aln.filter('c0i70 == id')

        if len(aln) > sample:
            keep = keep + aln.df[2:].id.sample(sample-2).tolist()
            aln = aln.filter(keep=keep) 

        if filename == None:
            filename = f'{model}.html'

        aln.to_html(consensus, filename, annotations=aln.df.iloc[0,0], remove_gaps=aln.df.iloc[1,0], fixed_index=True)

