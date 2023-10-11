# IPython log file

import os
import pathlib
import pandas as pd
from glob import glob
import rotifer
from rotifer.devel.beta import sequence as rdbs
logger = rotifer.logging.getLogger(__name__)
config = {
    "basedir": "/projects/salmonella/work/st/alndb",
    "pattern": "*.sto",
    "informat": "stockholm",
    "recursive": True,
    "ignore": [],
    "kwargs": {
        "pdb": {
            "basedir": "/projects/salmonella/work/st/colabfold",
            "pattern": "*_rank_001_*.pdb",
        },
    },
}

# Loading our alignments
class SequenceCollection():
    def __init__(
            self,
            basedir   = config["basedir"],
            pattern   = config["pattern"],
            informat  = config["informat"],
            recursive = config["recursive"],
            ignore    = config["ignore"],
            **kwargs
        ):
        if "kwargs" in config and config["kwargs"]:
            kwargs = { **config["kwargs"], **kwargs }
        if isinstance(pattern,str):
            pattern = [pattern]
        basedir = pathlib.Path(basedir)
        if not basedir.exists():
            logger.error(f'No directory named {basedir.name}')
            return None

        # Find alignments
        alndb = []
        alignments = dict()
        for extension in pattern:
            # Load data
            if recursive:
                it = basedir.rglob(extension)
            else:
                it = basedir.glob(extension)
            for x in it:
                if extension[0] == "*":
                    suffix = extension[1:]
                else:
                    suffix = extension
                basename = x.name.replace(f"{suffix}","")
                if basename in ignore:
                    continue
                logger.info(f'Loading {informat} file at {x.name}')
                alignments[basename] = rdbs.sequence(x.as_posix(),informat)

                # Load other paths
                otherfiles = []
                for colname, other in kwargs.items():
                    patt = os.path.join(basename, other["pattern"])
                    basepath = pathlib.Path(other["basedir"]) if "basedir" in other else basedir
                    if recursive:
                        files = list(basepath.rglob(patt))
                    else:
                        files = list(basepath.glob(patt))
                    files = [ y.as_posix() for y in files ]
                    if len(files) == 0:
                        files = None
                    elif len(files) == 1:
                        files = files[0]
                    otherfiles.append(files)

                # Add data to stack
                alndb.append((basename, x.as_posix(), *otherfiles))

        # Create internal dataframe for paths
        alndb = pd.DataFrame(alndb, columns=['name','path'] + list(kwargs.keys()))
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

    def view(self, model, *args, **kwargs):
        aln = self[model].copy()

        if "pdb_file" in kwargs:
            pdb = kwargs.pop("pdb_file")
        elif "pdb" in self.df.columns:
            pdb = self.df.loc[model,"pdb"]
            if isinstance(pdb,list):
                pdb = pdb[0]
        if pdb:
            aln = aln.add_pdb(aln.df.iloc[0,0], pdb_file=pdb)

        aln.view(*args, **kwargs)

    # Saving to HTML
    def to_html(self, model, filename=None, sample=50, consensus=90):
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

