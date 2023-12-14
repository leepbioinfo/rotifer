import os
import hashlib
import pathlib
import pandas as pd
import rotifer
logger = rotifer.logging.getLogger(__name__)
_defaults = {
    "basedir": "/databases/fadb",
    "checksum": True,
    "pattern": "*.fa",
    "recursive": True,
    "ignore": [],
}
config = loadConfig(__name__, defaults = _defaults)

# Classes
# Loading our alignments
class FileCollection():
    def __init__(
            self,
            basedir   = config["basedir"],
            pattern   = config["pattern"],
            recursive = config["recursive"],
            checksum  = config["checksum"],
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
        db = []
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
                db.append((basename, x.as_posix(), *otherfiles))

        # Create internal dataframe for paths
        colnames = ['name','path']
        if checksum:
            colnames.append('checksum')
        db = pd.DataFrame(db, columns=colnames + list(kwargs.keys()))

        # Store data in object
        self.df = db

    def __setitem__(self,key,value):
        self.df.loc[key] = value

    def __getitem__(self,key):
        return self.df.loc[key]

    def keys(self):
        return self.df.index.tolist()

# Is this library being used as a script?
if __name__ == '__main__':
    pass
