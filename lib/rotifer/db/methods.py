import types
import pandas as pd
from copy import deepcopy
import rotifer
logger = rotifer.logging.getLogger(__name__)

class GenomeCursor:
    def getids(self,obj):
        if isinstance(obj,types.NoneType):
            return set()
        elif isinstance(obj,set):
            return deepcopy(obj)
        elif not isinstance(obj,list):
            obj = [obj]
        assemblies = set()
        for s in obj:
            if hasattr(s,"assembly"):
                if isinstance(s.assembly,str):
                    assemblies.add(s.assembly)
                else:
                    logger.warn(f'Unknown assembly type {type(assembly)}: {assembly}')
            elif hasattr(s,"dbxrefs") and isinstance(s.dbxrefs,list):
                for x in s.dbxrefs:
                    if 'Assembly:' in x:
                        assemblyID = x.split(':')[-1]
                        assemblies.add(assemblyID)
        return assemblies

    def fetcher(self, accession):
        tries = self.tries
        targets = self.parse_ids(accession)
        stream = []
        for acc in targets:
            fh = self.open_genome(acc)
            if fh == None:
                continue
            fh.assembly = acc
            stream.append(fh)
        self.tries = tries
        return stream

    def parser(self, stream, accession):
        from Bio import SeqIO
        stack = []
        for fh in stream:
            if fh == None:
                continue
            for s in SeqIO.parse(fh,"genbank"):
                s.assembly = fh.assembly
                stack.append(s)
            fh.close()
        return stack

class GenomeFeaturesCursor(GenomeCursor):
    def getids(self,obj):
        if isinstance(obj,types.NoneType):
            return set()
        elif isinstance(obj,set):
            return deepcopy(obj)
        elif isinstance(obj,list):
            return set([ x.assembly for x in obj ])
        elif isinstance(obj,pd.DataFrame) and "assembly" in obj.columns:
            return set(obj.assembly)
        else:
            raise TypeError(f'Unknown object type {type(obj)}: {obj}')

    def parser(self, stream, accession):
        from Bio import SeqIO
        from rotifer.genome.utils import seqrecords_to_dataframe
        if not isinstance(stream, list):
            stream = [stream]
        data = []
        for fh in stream:
            datum = SeqIO.parse(fh,"genbank")
            datum = seqrecords_to_dataframe(
                datum,
                exclude_type = self.exclude_type,
                autopid = self.autopid,
                assembly = fh.assembly,
                codontable = self.codontable,
            )
            data.append(datum)
            fh.close()
        if len(data) > 0:
            data = pd.concat(data)
        else:
            data = seqrecords_to_dataframe([])
        return data

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch all accessions as a single dataframe.

        Parameters
        ----------
        accessions: list of strings
          Database accessions

        Returns
        -------
        Pandas dataframe or derived class.
        """
        from rotifer.genome.utils import seqrecords_to_dataframe
        stack = []
        for df in self.fetchone(accessions):
            stack.append(df)
        if stack:
            return pd.concat(stack, ignore_index=True)
        else:
            return seqrecords_to_dataframe([])
