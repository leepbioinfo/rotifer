import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import os
import re
logger = rotifer.logging.getLogger(__name__)

# Defaults
config = loadConfig(
    __name__.replace('rotifer.',':'),
    defaults = {
        'fetch': 'pfetch',
        'pdb_dir': os.path.join(os.environ['ROTIFER_DATA'] if 'ROTIFER_DATA' in os.environ else '/databases',"pdb"),
        'databases': ['pdb70','pfam'],
        'databases_path': os.path.join(os.environ['ROTIFER_DATA'] if 'ROTIFER_DATA' in os.environ else '/databases',"hhsuite"),
        'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
    })

class sequence:
    """
    This class represents a multiple sequence alignment (MSA)
    or a collection of unaligned sequences.

    It provides methods for MSA slicing, visualization, alignment
    rebuilding and annotation.

    Alignment annotations
    ---------------------
    Four types of annotations are supported by this class:

    - Row annotations
      Refers to properties of alignment rows, such as cluster
      assignments, functional annotations for entire sequences,
      such as EC numbers, or organism identifiers.

      Sequence annotations are stored as columns in the internal
      Pandas DataFrame and the user may manipulate such columns
      directly via the ``df`` attribute, as in the example below,
      where we remove the ``c80i30`` added by ``add_cluster``:

      >>> aln = aln.add_cluster(identity = 0.3)
      >>> aln.df.drop('c80i30', inplace=True, axis=1)

    - Character-based residue annotations
      Refers to a single column and are specified as strings of the
      same length as the alignment itself. Examples include
      **consensus sequences** and per-residue **secondary structure**
      codes.

      This type of annotation can be added to the internal Pandas
      DataFrame directly, i.e. using Pandas's ``append`` or
      ``insert`` methods but the user must be **careful** to set the
      value or the ``type`` column to something other than
      ``sequence``, which is researved for the aligned sequences.

    - Numerical values assigned to individual residues
      This type of annotation also assigns annotations to single
      columns of the alignment but are best described by numerical
      values, instead of single characters. Examples include the
      frequency of a given residue at each alignment column, the
      number of gaps per column or the column entropies.

    - Alignment or sequence regions
      Describe properties for subsets of columns in the alignment
      or for regions within a specific sequence. The former uses
      column indices to describe alignment features, while the later
      specifies a region within a given sequence using residue 
      coordinates that ignore gaps.

    Specialized methods are provided for the most common annotations
    (see, for example, ``add_cluster``, ``add_consensus``).

    Parameters
    ----------
    input_data : str or file
        Input file name, string or open file handler
    input_format : str
        Input alignment format. See Bio.SeqIO and/or
        Bio.AlignIO for a list of supported formats
    name:
        An (optional) name for the multiple sequence alignment

    See also
    --------
    Bio.SeqIO   : BioPython parser for aligned/unaligned sequences
    Bio.AlignIO : BioPython parser for aligned sequences
    pandas      : library used internally to store sequences

    Examples
    --------
    Simplest usage: visualize a (multi)FASTA file.

    >>> from rotifer.devel.beta.sequence import sequence
    >>> aln = sequence("alignment.aln")
    >>> aln.view()

    Load alignment from another file format, such as
    StockHolm and show the consensus above the alignment:

    >>> from rotifer.devel.beta.sequence import sequence
    >>> sto = sequence("alignment.sto", format="stockholm")
    >>> sto.add_consensus().view()

    Parsing a FASTA alignment from an explicit string:

    >>> from rotifer.devel.beta.sequence import sequence
    >>> aln = sequence(">Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n")
    >>> aln.add_consensus().view()
    """
    def __init__(self, input_data=None, input_format='fasta', name=None, local_database_path=config['local_database_path']):
        from io import IOBase
        self._reserved_columns = ['id','sequence','length','type']
        self.input_format = input_format
        self.name = name

        # Generate empty object
        if input_data is None:
            self.df = pd.DataFrame({}, columns=self._reserved_columns)
            self.input_format = None
            self.file_path = None

        # Initialize for strings
        elif isinstance(input_data,str):
            if os.path.exists(input_data):
                self.file_path = input_data
                if not name:
                    self.name = os.path.basename(input_data)
            else:
                self.file_path = 'StringIO'
                input_data = StringIO(input_data)
                if not name:
                    self.name = str(input_data)
            input_data = SeqIO.parse(input_data, input_format)
            self.df = self._seqrecords_to_dataframe(input_data)

        # Initialize for IO streams
        elif isinstance(input_data, IOBase):
            input_data = SeqIO.parse(input_data, input_format)
            self.df = self._seqrecords_to_dataframe(input_data)

        # Initialize for list of Bio.SeqRecords
        elif isinstance(input_data, list):
            if isinstance(input_data[0],SeqRecord):
                self.df = self._seqrecords_to_dataframe(input_data)
            elif isinstance(input_data[0],str):
                import rotifer.db as rdb
                seqrecords = rdb.proteins(input_data, local_database_path=local_database_path)
                self.df = self._seqrecords_to_dataframe(seqrecords)

        # Initialize for Pandas DataFrame
        elif isinstance(input_data, pd.DataFrame):
            other = [ x for x in input_data.columns if x not in self._reserved_columns ]
            self.df = input_data[['id','sequence']]
            self.df['type'] = 'sequence'
            self.df[other] = input_data[other]

        # Initialize for Pandas Series
        elif isinstance(input_data, pd.Series):
            self.df = input_data.dropna().reset_index()
            self.df.columns = ['id','sequence']
            self.df['type'] = 'sequence'

        # For creating an empty object
        else:
            self.df = pd.DataFrame(columns=self._reserved_columns)

        # Still unnamed? Try using the input filename
        if not name:
            if hasattr(input_data,"stream"):
                if hasattr(input_data.stream,"filename") and callable(input_data.stream.filename):
                    self.name = os.path.basename(input_data.stream.filename())
                elif hasattr(input_data.stream,"name"):
                    self.name = os.path.basename(input_data.stream.name)
            elif not self.df.empty:
                self.name = self.df.id.loc[0]

        # Make sure the new object is clean!
        self._reset()

    # Functionis to alow item assignment to the sequence class.
    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getitem__(self, key):
        return getattr(self, key)

    def _seqrecords_to_dataframe(self, data):
        cols = self._reserved_columns + ['description']
        parsed = []
        for s in data:
            current = [ s.id, str(s.seq) ]
            current.append(len(current[1].replace("-","")))
            current.append('sequence')
            current.append(s.description)
            parsed.append(current)
        df = pd.DataFrame(parsed, columns=cols)
        return df

    def _reset(self):
        # Recalculate each sequence's length
        maxlen = self.df.sequence.str.len().max()
        self.df['length'] = np.where(
                self.df.type == "sequence",
                self.df.sequence.str.replace('-', '').str.len(),
                maxlen)
        if not self.name and not self.df.empty:
            self.name = self.df.id.loc[0]
        self.df.reset_index(drop=True, inplace=True)

        # Statistics holder
        if not hasattr(self,'numerical'):
            self.numerical = pd.DataFrame(columns=['type']+list(range(1,self.get_alignment_length()+1)))

    def __len__(self):
        return len(self.df.query('type == "sequence"'))

    def copy(self, reset=False):
        '''
        Copy alignment to new object.
        '''
        from copy import deepcopy
        result = deepcopy(self)
        if reset:
            result._reset()
        return result

    def filter(self, query=None, minlen=0, maxlen=0, keep=[], remove=[], regex=None):
        '''
        Function to filter / select rows from the MSA.

        Parameters
        ----------
        query: str
          String to use as input for the ``query`` method of the
          internal Pandas Dataframe
        minlen : int
          Length of the shortest sequence, without gaps
        maxlen : int
          Length of the longest sequence, without gaps
        keep : list of strings
          List of sequence IDs to preserve

          The ``keep`` parameter has the highest precedence over
          any other parameters, i.e. any sequence listed in
          ``keep``, is **never** removed.

        remove : list of strings
          List of sequences IDs to be removed.

          Sequences in this list are removed even if failing the
          criteria set by other parameters, except ``keep``.

        regex : regular expression
          Regular expression to identify rows based on
          residue patterns matching their aminoacid/nucleotide
          sequences.

        Returns
        -------
        New MSA

        Examples
        --------
        Filter sequences with at least 50 and at most 100 residues

        >>> aln.filter(minlen=50, maxlen=100)

        Filter sequences based on a cluster but make sure NPE27555.1
        is also kept

        >>> aln = aln.add_cluster(identity=0.6)
        >>> aln.filter('c80i80 == "WP_091936315.1"', keep="NPE27555.1")
        '''
        result = self.copy(reset=True)

        # Build query statement
        querystr = []
        if query:
            querystr.append(query)
        if minlen > 0:
            querystr.append('(length >= @minlen)')
        if maxlen > 0:
            querystr.append('(length <= @maxlen)')
        if regex:
            querystr.append('(sequence.str.replace("-", "").str.contains(@regex, regex=True))')
        if remove:
            querystr.append('(id not in @remove)')
        querystr = " and ".join(querystr)
        if keep:
            if querystr:
                querystr = f'(id in @keep) or ({querystr})'
            else:
                querystr = f'id in @keep'
        if querystr:
            result.df = result.df.query(querystr)
        return result

    def get_alignment_length(self):
        '''
        Retrieve the number of columns in the alignment.

        Returns
        -------
        Integer

        Examples
        --------
        >>> al = aln.get_alignment_length()
        '''
        if self.df.empty:
            return 0
        else:
            return int(self.df.query('type == "sequence"').sequence.str.len().max())

    def slice(self, position):
        '''
        Select and concatenate one or more sets of columns.

        Coordinate systems
        ------------------
        This method cuts the input aligment based on two types of
        coodinate systems:

        * Alignment-based coordinates
          Refer to the column indices in the original alignment

        * Sequence-based coordinates
          Refer to residues in a reference sequence, without gaps

        All types of coordinates systems above are **one-based**
        (residue-based) and closed on both ends, i.e. the first
        column or residue position is 1 and the last column equals
        the length of the alignment or sequence.

        Notice that coordinates must be provided as tuples (see below).

        Examples
        --------

        - Alignment-based coordinates
          Fetch columns 20 to 130:

          >>> aln.slice((20,130))

        - Sequence-based coordinates
          Fetch columns corresponding to residues 10 to 80 of
          the sequence WP_003247817.1

          >>> aln.slice((10,80,"WP_003247817.1"))

        - You can also fetch several non-overlapping regions
          and these will be concatenated after removal of the
          interviening regions:

          >>> aln.slice([ (10,80),(110,160) ])

        '''
        result = self.copy()
        if isinstance(position[0],int):
            position = [position]

        # Cut slices
        sequence  = []
        numerical = []
        for pos in position:
            pos = [*pos]
            if len(pos) == 3:
                refseq = pd.Series(list(result.df.query(f'''id == "{pos[2]}"''').sequence.values[0]))
                refseq = refseq.where(lambda x: x != '-').dropna().reset_index().rename({'index':'mapped_position'}, axis=1)
                pos[0:2] = (refseq.loc[pos[0]-1:pos[1]].mapped_position.agg(['min','max'])).tolist()
                pos[0] += 1
            sequence.append(result.df.sequence.str.slice(pos[0]-1, pos[1]))
            numerical.extend(list(range(pos[0],pos[1]+1)))

        # Rebuild sequence
        result.df['sequence'] = pd.concat(sequence, axis=1).sum(axis=1)

        # Rebuild numerical
        other = set(self.numerical.columns) - set(['type']) - set(list(range(1, self.get_alignment_length() + 1)))
        result.numerical = result.numerical[['type'] + numerical + list(other)]
        result.numerical.columns = ['type'] + list(range(1,len(numerical)+1)) + list(other)

        # Return new sequence object
        result._reset()
        return result

    def sort(self, by=['length'], ascending=True, inplace=False, id_list=None, tree_file=None, tree_format='newick'):
        """
        Sort alignment rows.

        Parameters
        ----------
        by : list of strings, default ['length']
            List of criteria for sorting.

            The following criteria are always supported:
            - id : sort by sequence identifiers
            - length : sequence length
            - list : order follows a list of sequence IDs
                See parameter ``id_list``.
            - name : sort by sequence identifiers
            - tree : sort by phylogeny
                See parameter ``tree_file``.
                Leaf names must match sequence identifiers.
            - residues: alignment columns

            Additionally, user-supplied sequence annotations may
            also be used to sort alignment rows. See the section
            ```Alignment annotations```.

        ascending : bool or list of bools, default True
            If True, sort in ascending order, otherwise descending.

        inplace : bool, default False
            Sort in place, i.e. without copying of the object

        id_list : list of strings or file name
            A list of sequence identifiers to use in sorting.
            If a file name is given as a string, the list will be
            automatically loaded from that file.

        tree_file : path to phylogeny
            For this parameter to be effective, at least a subset
            of the leaves in the phylogeny must match sequence
            identifiers in the alignment.

            Leaf names that do not match sequence identifiers are
            ignored. Sequence identifiers from the alignment that
            are absent in the phylogeny are listed after

        tree_format: phylogeny file format
            All formats supported by Ete3 or BioPython are valid

        Returns
        -------
        MSA with rows in a new order

        See also
        --------
        sequence.annotations : 

        Examples
        --------
        Default: sort by sequence length, ignoring gaps,
        from shortest to longest:

        >>> aln2 = aln.sort()

        Sort by length and then by the sequence cluster
        (see ``add_cluster``).

        >>> aln = aln.add_cluster(coverage=0.8, identity=0.3)
        >>> aln = aln.sort(['length','c80i30'])

        Make sure sequence ``PUA33204.1`` shows up first.
        Sort the others by length.

        >>> aln.sort(['list','length'], id_list=['PUA33204.1']).view()

        Sort in-place by columns 32 and 27 and by length:

        >>> aln.sort([32,27,'length'], inplace=True)

        """
        import sys

        # Inplace mode
        if inplace:
            result = self
        else:
            result = self.copy()
        result.df.reset_index(inplace=True, drop=True)

        # Prepare local variables
        cols = result.df.columns
        isnumeric = [ x for x in by if x not in cols and isinstance(x,int) ]
        if isnumeric:
            matrix = self.residues

        fields = []
        identifiers = result.df.id
        for item in by:
            if item == 'tree':
                from ete3 import Tree
                tree = Tree(tree_file)
                leaves = pd.Series([ x.name.strip("'") for x in tree.get_leaves() ], name='leaves')
                leaves = leaves.reset_index().rename({'index':'order'}, axis=1).set_index('leaves').order.to_dict()
                fields.append(identifiers.where(lambda x: x.isin(leaves)).replace(leaves).rename("tree"))

            elif item == 'list':
                if isinstance(id_list,str) and os.path.exists(id_list):
                    ids = pd.read_csv(id_list, names=['id']).id
                elif isinstance(id_list, list):
                    ids = pd.Series(id_list, name='id')
                ids = ids.reset_index().rename({'index':'order'}, axis=1).set_index('id').order.to_dict()
                fields.append(identifiers.where(lambda x: x.isin(ids)).replace(ids).rename("list"))

            elif item == 'name':
                fields.append(identifiers.rename("name"))

            elif item in cols:
                fields.append(result.df[item])

            elif isinstance(item,int):
                fields.append(matrix.loc[:,item])

            else:
                logger.error(f'sequence.sort: Unsupported criteria or missing annotation ({item})')

        # Concatenate sorting fields
        fields = pd.concat(fields, axis=1).sort_values(by=by, ascending=ascending)
        #result.df = result.df.reindex(fields.index)
        result.df = result.df.loc[fields.index]
        if not inplace:
            return result

    def add_cluster(self, coverage=0.8, identity=0.7, name=None, inplace=False):
        '''
        Add or update MMseqs2 clustering data for all sequences.

        Parameters
        ----------
        coverage : float
            Minimum coverage required for both sequences in each
            pairwise alignment.
        identity : float
            Minimum percentage identity required for both sequences
            in each pairwise alignment.
        name: string
            Set the name of the annotation column for the clusters
        inplace: bool
            Add column inplace, without creating a copy of the MSA
            If set to True, this method return nothing.

        Returns
        -------
        A new MSA object

        Examples
        --------
        Add an annotation column listing representatives of clusters
        defined by 60% alignment coverage and 30% identity for both
        sequences

        >>> aln.add_cluster(coverage=0.6, identity=0.3)
        '''
        import tempfile
        from subprocess import Popen, PIPE, STDOUT

        # Adjust parameters
        if identity > 1:
            identity /= 100
        if coverage > 1:
            coverage /= 100
        if not name:
            name = f'c{int(float(coverage)*100)}i{int(float(identity)*100)}'

        path = os.getcwd()
        result = None
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            self.to_file(f'{tmpdirname}/seqaln', remove_gaps=True)
            Popen(f'mmseqs easy-cluster {tmpdirname}/seqaln nr --min-seq-id {identity} -c {coverage} tmp', stdout=PIPE,shell=True).communicate()
            d = pd.read_csv(f'{tmpdirname}/nr_cluster.tsv', sep="\t", names=['cluster', 'pid']).set_index('pid').cluster.to_dict()
            if inplace:
                self.df[name] = self.df.id.replace(d)
            else:
                result = self.copy()
                result.df[name] = result.df.id.replace(d)
        os.chdir(path)

        return result

    def add_consensus(self, cutoffs=(50, 60, 70, 80, 90), separator=None):
        """
        Add consensus rows to the MSA object.

        Parameters
        ----------
            cutoffs : list of integers
                Minimum frequency, in each column, required for
                the conserved residue or residue category
            separator : character
                Add a row with filled with the user-defined
                character to separate the alignment and the
                consensus lines

        Returns
        -------
        A copy of the MSA with consensus sequences
        """
        result = self.copy()
        cx = []
        frequencies = result.residue_frequencies
        for cutoff in reversed(sorted(cutoffs)):
            consensus = self.consensus(cutoff, frequencies=frequencies)
            cx.append([ f'Consensus {cutoff}%', consensus, len(consensus.replace('.','')), "consensus" ])
        if separator:
            cx.append([ 'separator', "".join([ separator for x in range(0,self.get_alignment_length()) ]), self.get_alignment_length(), "view" ])
        result.df = pd.concat([pd.DataFrame(cx, columns=self._reserved_columns), result.df])
        return result

    def add_hhpred(self, hhpred_data, hhpred_id, hhsearch=False):
        '''
        Merge the MSA and a pairwise alignment generated by HH-suite.

        The resulting MSA will include both the sequence and 
        secondary structure of the selected HH-suite result.

        Parameters
        ----------
        hhpred_data : string or file path
            HH-suite output file or string
        hhpred_id : string
            HH-suite hit identifier
        '''
        import re

        # Parser
        try:
            with open(hhpred_data, 'r') as f:
                texto = f.read()
        except:
            texto = hhpred_data

        # Parse
        match = re.findall(f'No {hhpred_id}(.+?)No ', texto, re.DOTALL)[0].strip()
        target = re.findall('T Consensus.*?\nT\s(.*?)\s',match, re.MULTILINE)[0]
        if hhsearch:
            query = re.findall('Query\s+(.*?)\s',texto, re.MULTILINE)
        else:
            query = re.findall('Q ss_pred.*?\nQ\s(.*?)\s',match, re.MULTILINE)
        
        if query:
            query = query[0]
        T_con = ''.join(re.findall('T Consensus.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        sequence_target = ''.join(re.findall(f'T\s+{target}.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        sequence_query = ''.join(re.findall(f'Q\s+{query}.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        T_ss_pred = ''.join(re.findall('T ss_pred\s+(.*?)$',match,  re.MULTILINE))
        T_ss_dssp = ''.join(re.findall('T ss_dssp\s+(.*?)$',match,  re.MULTILINE))
        Q_ss_pred = ''.join(re.findall('Q ss_pred\s+(.*?)$',match,  re.MULTILINE))
        Q_con = ''.join(re.findall('Q Consensus.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        if hhsearch:
            structure_df = pd.DataFrame({'query':list(sequence_query), 'target':list(sequence_target), 'ss':list(T_ss_dssp)})
        else:
            structure_df = pd.DataFrame({'query':list(sequence_query), 'query_pred': list(Q_ss_pred), 'target':list(sequence_target), 'ss':list(T_ss_dssp)})
        
        # join aln to hhpred results
        result = self.copy()
        aln = pd.Series(list(result.df.query('id == @query').sequence.iloc[0])).where(lambda x: x!='-').dropna().reset_index().rename({'index':'position', 0:'sequence'}, axis=1)
        s_aln = ''.join(aln.sequence).find(''.join(structure_df['query'].where(lambda x : x !='-').dropna()).upper())
        e_aln = s_aln + len(''.join(structure_df['query'].where(lambda x : x !='-').dropna())) -1
        aln.loc[s_aln : e_aln , 'dssp'] = structure_df[structure_df['query'] != '-'].ss.to_list()
        if not hhsearch:
            aln.loc[s_aln : e_aln , 'sspred'] = structure_df[structure_df['query'] != '-'].query_pred.to_list()
        aln.loc[s_aln : e_aln , target] = structure_df[structure_df['query'] != '-'].target.to_list()
        aln = pd.Series(list(result.df.iloc[0].sequence)).to_frame().join(aln.set_index('position', drop=True)).fillna('-')
        aln = aln.T.apply(lambda x: ''.join(list(x)), axis=1).reset_index().rename({'index': 'id', 0: 'sequence'}, axis=1) 
        result.df =  pd.concat([aln,result.df]).iloc[2:].reset_index(drop=True)

        return result

    def add_pdb(self, pdb_id, chain_id='A', pdb_file=None, pdb_dir=config['pdb_dir']):
        """
        Add structure annotation from PDB.

        Parameters
        ----------
        pdb_id : string
            4-letter PDB code or arbitrary string
        chain_id : string, defaul 'A'
            Target chain identifier.
        pdb_file : string
            PDB file or URL of remote PDB file.
        pdb_dir : string
            Path to a local PDB mirror or a directory
            where PDB files are stored

        See also
        --------
        Rotifer's configuration
        """
        import tempfile
        from Bio import pairwise2
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        import Bio.PDB.PDBList as PDBList
        result = self.copy()

        # Find local file or download and then open it!
        if pdb_file:
            if not os.path.exists(pdb_file):
                if os.path.exists(os.path.join(pdb_dir,pdb_file)):
                    pdb_file = os.path.exists(os.path.join(pdb_dir,pdb_file))
                else:
                    # Try using pdb_file as URL
                    import urllib
                    pdb_data = urllib.request.urlopen(pdb_file).read()
                    pdb_file = tempfile.NamedTemporaryFile(suffix=".pdb", delete=True)
                    pdb_file.write(pdb_data)
                    pdb_file.flush()
                    pdb_file.seek(0)
        else:
            # No file!
            if os.path.exists(os.path.join(pdb_dir,pdb_id[1:3].lower(),"pdb"+pdb_id.lower()+".ent.gz")):
                # Search PDB code in local PDB mirror
                pdb_file = os.path.join(pdb_dir, pdb_id[1:3].lower(), "pdb"+pdb_id.lower()+".ent.gz")
            else:
                pdb_file = PDBList(verbose=False).retrieve_pdb_file(pdb_id.lower(), file_format='pdb', pdir=GlobalConfig['cache'])
        if isinstance(pdb_file,str) and os.path.exists(pdb_file):
            if pdb_file[-3:] == '.gz':
                import gzip
                orig = gzip.open(pdb_file,"rt")
                pdb_file = tempfile.NamedTemporaryFile(suffix=".pdb", delete=True, mode="r+t")
                pdb_file.write(orig.read())
                pdb_file.flush()
                pdb_file.seek(0)
                orig.close()
            else:
                pdb_file = open(pdb_file,"rt")

        # Parse PDB file to a Pandas DataFrame
        dssp_columns  = ['pdb_id','chain','c1','c2','c3','idx','aa','secstr','rASA','phi','psi']
        dssp_columns += ['NH_O_1_relidx','NH_O_1_energy','O_NH_1_relidx','O_NH_1_energy']
        dssp_columns += ['NH_O_2_relidx','NH_O_2_energy','O_NH_2_relidx','O_NH_2_energy']
        dssp = PDBParser().get_structure(pdb_id, pdb_file)
        dssp = DSSP(dssp[0], pdb_file.name, file_type="PDB")
        dssp = pd.DataFrame([ (pdb_id,x[0],*x[1],*dssp[x]) for x in dssp.keys() ], columns=dssp_columns)
        dssp_to_ehc = {"I":"C","S":"C","H":"H","E":"E","G":"C","B":"E","T":"C","C":"C"}
        dssp.secstr = dssp.secstr.map(dssp_to_ehc).fillna('-')
        dssp = dssp.query('chain == @chain_id')
        pdb_file.close()
        if os.path.exists(pdb_file.name) and pdb_file.name[0:len(GlobalConfig['cache'])] == GlobalConfig['cache']:
            os.remove(pdb_file.name)

        # Search for pdb sequence in the aligment
        pdb_from_pdb = ''.join(dssp.aa.to_list())
        pdb_index = result.df[result.df.id.str.contains(pdb_id, case=False)]
        pdb_index = pdb_index[~pdb_index.id.str.contains('ss_from')] #Avoids annotating itself
        pdb_index = pdb_index.index[0]
        pdb_in_aln = result.df.loc[pdb_index].sequence.replace('-', '')
        ali = pairwise2.align.localxx(pdb_in_aln, pdb_from_pdb)[0]
        ss = pd.Series(list(ali.seqB)).where(lambda x : x != '-').dropna().to_frame()
        ss['structure'] = dssp.secstr.to_list()
        pdb_df = pd.Series(list(ali.seqA)).rename('seq').to_frame().join(ss['structure']).fillna('-').query(' seq != "-"')
        to = pd.Series(list(result.df.loc[pdb_index].sequence)).where(lambda x: x !='-').dropna().rename('ung').to_frame()
        to['structure'] = pdb_df.structure.to_list()
        pdbn = f'ss_from:{pdb_id}_{chain_id}'
        pdbss = ''.join(pd.Series(list(result.df.loc[pdb_index].sequence)).to_frame().join(to).fillna('-').structure.to_list())
        result.df =  pd.concat([pd.DataFrame([[pdbn,pdbss,len(pdbss),"residue_annotation"]], columns=self._reserved_columns),self.df])
        return result

    def add_seq(self, seq_to_add, cpu=12, fast=False, fetch=config["fetch"]):
        import tempfile
        import subprocess
        from rotifer.db.ncbi import NcbiConfig
        from subprocess import Popen, PIPE, STDOUT

        # Make sure input is a list
        if not isinstance(seq_to_add,list):
            seq_to_add = [ seq_to_add ]

        # Dump input sequences
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Save current sequences to temporary file
            self.to_file(f'{tmpdirname}/seqaln') 

            # Save new sequences to temporary file
            from Bio.SeqRecord import SeqRecord
            if isinstance(seq_to_add[0],SeqRecord):
                SeqIO.write(seq_to_add, f'{tmpdirname}/acc.fa', "fasta")
            else:
                pd.Series(seq_to_add).to_csv(f'{tmpdirname}/acc', index=None, header=None)
                Popen(f'{fetch} {tmpdirname}/acc > {tmpdirname}/acc.fa' , stdout=PIPE,shell=True).communicate()

            # Run MAFFT
            child = f'mafft --thread {cpu} --add {tmpdirname}/acc.fa'
            if not fast:
                child = f'{child} --maxiterate 1000 --localpair'
            child = Popen(f'{child} {tmpdirname}/seqaln', stdout=PIPE, shell=True).communicate()

        # Load output alignment
        result = self.from_string(child[0].decode("utf-8"), input_format = 'fasta')
        result.file_path = 'from add_seq function'
        return result

    def consensus(self, cutoff=50, frequencies=None):
        '''
        Generate consensus strings for the alignment object

        Parameters
        ----------
        cutoff : integer or float, default 50
          Minimum frequency for conserved residues or categories
        frequencies: Pandas DataFrame
          Frequencies of residues at each column
          See sequence.residue_frequencies

        Returns
        -------
        A string representing the consensus sequence
        '''

        # Ranking of amino acid categories
        aa_type_dict =  {'a': 6,
            'l':4,
            'h':8,
            '+':3,
            '-':1,
            'c':7,
            'p':10,
            'o':0,
            'u':5,
            's':9,
            'b':11,
            '_':12,
            '.':13,
            'gap':14
        }

        # Copying frequency table and building consensus
        if isinstance(frequencies,pd.DataFrame) and not frequencies.empty:
            result = frequencies
        else:
            result = self.residue_frequencies
        result.rename({'gap':'.'}, inplace=True)
        result = pd.concat([result, pd.DataFrame(columns=result.columns, index=['.']).fillna(cutoff+100)])
        result = result.melt(ignore_index=False).reset_index().rename({'index':'aa', 'variable':'position', 'value':'freq'}, axis=1)
        result['ranking'] = result.aa.map(aa_type_dict)
        result = result.query(f'freq >= {cutoff}').sort_values(['position','ranking'], na_position='first').drop_duplicates(subset='position')
        return ''.join(result.aa.to_list())

    @property
    def freq_table(self):
        '''
        DEPRECATED: use residue_frequencies instead!
        Residue and residue categories frequencies.
        '''
        return self.residue_frequencies

    @property
    def residues(self):
        '''
        Access alignment columns as a Pandas dataframe.
        '''
        return self.df.query('type == "sequence"').sequence.str.split("", expand=True).loc[:,1:self.get_alignment_length()]

    @property
    def residue_frequencies(self):
        '''
        Column-wise frequencies of residues and residue categories.
        '''
        freq_df = self.residues
        freq_df = freq_df.apply(pd.value_counts).fillna(0).astype(int)/len(freq_df)*100 
        freq_df.rename({'-':'gap','.':'gap','?':'X'}, inplace=True)

        aromatic = ['F','Y', 'W', 'H']
        alifatic = ['I','V','L']
        hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
        positive = ['H', 'K', 'R']
        negative = [ 'D', 'E']
        charged = positive + negative
        polar = charged + ['Q', 'N', 'S', 'T','C']
        alcohol = ['S','T']
        tiny = ['G', 'A', 'S']
        small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
        big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
        all_aa = ['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T']

        aa_type_names =  {'aromatic':[aromatic,'a'],
                         'alifatic':[alifatic,'l'],
                         'hydrophobic': [hydrophobic,'h'],
                         'positive':[positive,'+'], 
                         'negative':[negative,'-'],
                         'charged':[charged,'c'], 
                         'polar':[polar,'p'],
                         'alcohol':[alcohol,'o'],
                         'tiny':[tiny,'u'],
                         'small':[small,'s'],
                         'big':[big,'b'],
                         'all_aa':[all_aa, '_']
                         }

        missing = set(all_aa) - set(freq_df.index.tolist())
        for x in missing:
            freq_df.loc[x] = 0
        for x in aa_type_names.keys():
            freq_df = pd.concat([freq_df,pd.DataFrame({aa_type_names[x][1]:freq_df.loc[aa_type_names[x][0]].sum()}).T])
        return freq_df

    def to_a3m(self, name=None, remove_gaps=False):
        name = self.name
        if not name:
            name = self.df.id.iloc[0]
        a3m = f'#{name}\n'
        ss = self.df.query('type == "structure prediction"')
        if not ss.empty:
            a3m += ">ss_pred " + ss.query('id != "conf"').id.iloc[0].upper() + "\n"
            a3m += ss.query('id != "conf"').sequence.iloc[0].replace(" ","-") + "\n"
            a3m += ">ss_conf " + ss.query('id != "conf"').id.iloc[0].upper() + "\n"
            a3m += ss.query('id == "conf"').sequence.iloc[0].replace(" ","0") + "\n"
        a3m += self.to_string(remove_gaps=remove_gaps)
        return a3m

    def to_bioalign(self):
        """
        Convert the MSA to BioPython's Bio.Align.MultipleSeqAlignment.

        See also
        --------
        Bio.Align.MultipleSeqAlignment : MSA object in BioPython
        Bio.AlignIO : BioPython I/O modules for alignments
        """
        from Bio.Align import MultipleSeqAlignment
        return MultipleSeqAlignment(self.to_seqrecords())

    def to_color(self, color='fg', scale=True, interval=10):
        """
        Colorize sequence column and return aligment as a DataFrame
        """

        def color_res(s, cs):
            if s in 'A I L M F W V'.split():
                return cs(s, "033")
            elif s in 'K R'.split():
                return cs(s, "124")
            elif s in 'E D'.split():
                return cs(s, "127")
            elif s in 'N Q S T'.split():
                return cs(s, "034")
            elif s  == 'C':
                return cs(s, "168")
            elif s == 'G':
                return cs(s, "166")
            elif s == 'P':
                return cs(s, "178")
            elif s in 'H Y'.split():
                return cs(s, "037")
            else:
                 return cs(s, "000")

        def color_bg(s, color = ''):
            '''
            s: String
            '''
            if color == "000":
                color = f'49;5;000'
            else:
                color = f'48;5;{color}'
            return f'\033[{color}m{s}\033[0m'

        def color_fg(s, color = ''):
            if color == "000":
                color = f'39;5;000'
            else:
                color = f'38;5;{color}'
            return f'\033[{color}m{s}\033[m'

        # Make a copy of the input dataframe
        result = self.df.copy()
        header = pd.Series(result.columns, index=result.columns).to_frame().T

        # Add scale
        if scale:
            scale = self._scale_bar(self.get_alignment_length(), interval=interval)
            result = pd.concat([ header, scale, result.astype(str) ])
        else:
            result = pd.concat([ header, result.astype(str) ])

        # Color the sequence column and return
        color_switch = {'background':color_bg, 'bg':color_bg, 'foreground':color_fg, 'fg':color_fg}
        result.sequence = result.sequence.str.pad(result.sequence.str.len().max(), side="right")
        result.sequence = result.sequence.map(lambda x: ''.join([color_res(y, color_switch[color]) for y in x]))
        return result

    def to_file(self, file_path=None, output_format='fasta', annotations=None, remove_gaps=False):
        if output_format == "a3m":
            fh = open(file_path,"wt")
            fh.write(self.to_a3m())
            fh.flush()
            fh.close()
            return len(self)
        else:
            return SeqIO.write(self.to_seqrecords(annotations=annotations, remove_gaps=remove_gaps), file_path, output_format)

    def to_hist(self, bins=10):
        """
        This method displays the distribution of sequence
        lengths for all the sequences in the alignment, 
        after removal of all their gaps.

        Parameters
        ----------
        bins : int
            Number of histogram bins to discretize the distribution

        Returns
        -------
        None

        Examples
        --------
        >>> aln.hist(50)
        """
        from ascii_graph import Pyasciigraph
        import collections
        collections.Iterable = collections.abc.Iterable
        print(f'Total proteins: {len(self.df)}')
        a = self.df['length'].value_counts().to_frame().reset_index()
        a = a.sort_values('index')
        a['raw_bin'] = pd.cut(a['index'],bins,precision=0)
        a['bin'] = a.raw_bin.apply(lambda x : '{} - {}'.format(int(x.left),int(x.right)))
        test = a.groupby('bin').agg({'length':'sum'}).reset_index().apply(tuple, axis=1)
        graph = Pyasciigraph()
        for line in  graph.graph('count \t sequence size', test):
            print(line)

    def to_seqrecords(self, annotations=None, remove_gaps=False):
        """
        Convert the MSA to a list of BioPython's SeqRecord objects.
        """
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # Access and filter internal dataframe
        df = self.df
        if annotations:
            df = df.query('type == "sequence" or type in @annotations').reset_index(drop=True)
        else:
            df = df.query('type == "sequence"')

        # Convert each row to a SeqRecord
        result = []
        for row in list(df.T.to_dict().values()):
            if remove_gaps:
                row["sequence"] = row["sequence"].replace("-","")
            if "description" not in row:
                row["description"] = "unknown sequence"
            if not isinstance(row["description"], str):
                row["description"] = "unknown sequence"
            result.append(SeqRecord(id=row["id"], seq=Seq(row["sequence"]), description=row["description"]))
        return result

    def to_string(self, output_format='fasta', annotations=None, remove_gaps=False):
        sio = StringIO("")
        SeqIO.write(self.to_seqrecords(annotations=annotations, remove_gaps=remove_gaps), sio, output_format)
        return sio.getvalue()

    def align(self, method='famsa', cpu=12, region=False):
        """
        Rebuild the alignment using Mafft.

        Parameters
        ----------
        method : string, default is famsa
            famsa ,
            mafft runs maaft with defaul parameter 
            linsi runs mafft  using the Smith-Waterman algorithm with options ```--maxiterate 1000``` and ```--localpair```
            kalign 
        cpu : integer, default is 12
            Number of threads to use, kalign do not have thread option

        Returns
        -------
        A new MSA object.
        """
        from subprocess import Popen, PIPE, STDOUT
        seq_string = self.to_string(remove_gaps=True).encode()
        if region:
            seq_string = self.slice((region[0],region[1])).to_string(remove_gaps=True).encode()

        if method =='mafft':
            child = Popen(f'cat|mafft  --thread {cpu} -' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        elif method == 'linsi':
            child = Popen(f'cat|mafft  --maxiterate 1000 --localpair --thread {cpu} -' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        elif method =='famsa':
            child = Popen(f'cat|famsa -t {cpu} -v STDIN STDOUT' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        elif method =='kalign':
            child = Popen(f'cat|kalign' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)

        result = self.from_string(child[0].decode("utf-8"), input_format = 'fasta')
        if region:
            r1 = self.slice((1, region[0] - 1))
            r2 = self.slice((region[1] + 1, len(self.df.iloc[0,1])))
            r3 = self.copy()
            r3.df.sequence = r1.df.sequence + result.df.sequence + r3.df.sequence
            result = r3

        result.file_path = 'from align function'
        result.name = self.name
        return result

    def view(self, color=True, scale=True, consensus=True, separator="=", interval=10, columns=True, pager='less -SR'):
        """
        Display alignment and alignment annotations.

        Parameters
        ----------
        color : bool
            Color sequence residues
        scale : bool, default True
            If set to False, no scale is shown
        consensus : bool, default True
            Whether to display consensus rows
        separator : character
            Single character to fill row separating the alignment
            and the consensus sequence
        interval : integer, default 10
            Interval between position marks in the scale
        columns : bool or list of strings, default is True
            List of annotation columns to show
            If set to False, only the default columns are shown
            If set to True, all columns in the internal DataFrame are shown
        pager: string
          Command line to be used for displaying the alignment.

        Returns
        -------
            None
        """
        from IPython.core.page import page
        df = self.copy()
        basic_columns = ['id', 'sequence','length','type',]
        if isinstance(columns,bool):
            if not columns:
                df.df = df.df[basic_columns]
        else:
            if isinstance(columns,list):
                df.df = df.df[basic_columns + columns]
            else :
                logger.error('columns should be either a list or a bool')

        if consensus:
            df = df.add_consensus(separator=separator)
        if color:
            df = df.to_color(scale=scale, interval=interval)
            page(df.to_string(header=False, index=False) + "\n", pager_cmd=pager)
        else:
            df = df.df.copy()
            if scale:
                scale = self._scale_bar(self.get_alignment_length(), interval=interval)
                df = pd.concat([ scale, df ]) 
                df.sequence = df.sequence.str.pad(df.sequence.str.len().max(), side="right")
            page(df.to_string(index=False) + "\n", pager_cmd=pager)

    def hhblits(self, databases=config['databases'], database_path=config['databases_path'], view=True):
        """
        Search the alignment against a HMM databases using hhsearch.

        Parameters
        ----------
        databases : list of strings, default ['pfam','pdb70']
            List of HMM databases to include in the search
        database_path : string, default is ROTIFER_DATA/hhsuite
            Path to the directory where the HMM databases are stored
        view : bool, default True

        Returns
        -------
            A tuple of two elements:
            - The HHblits output as a string
            - HHblits output as a Pandas DataFrame
              See rotifer.io.hhsuite

        See also
        --------
            rotifer environment configuration

        Examples
        --------
        Load alignment in multi-FASTA format and compare it to Pfam

        >>> aln = sequence("myaln.aln")
        >>> (hhout,hhdf) = aln.hhblits(databases=['pfam'])
        """
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer.io.hhsuite import read_hhr

        dbs = " ".join([ " -d " + os.path.join(database_path, x) for x in databases ])
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/seqaln')
            child = f'hhblits -i {tmpdirname}/seqaln {dbs} -M 50 -cpu 18 -o {tmpdirname}/seqaln.hhr'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            hhtable = read_hhr(f'{tmpdirname}/seqaln.hhr')
            with open(f'{tmpdirname}/seqaln.hhr') as f:
                hhr_result = f.read()
            if view:
                from IPython.core.page import page
                page(hhr_result)
            return (hhr_result, hhtable)

    def hhsearch(self, databases=config['databases'], database_path=config['databases_path'], view=True):
        """
        Search the alignment against a HMM databases using hhsearch.

        Parameters
        ----------
        databases : list of strings, default ['pfam','pdb70']
            List of HMM databases to include in the search
        database_path : string, default is ROTIFER_DATA/hhsuite
            Path to the directory where the HMM databases are stored
        view : bool, default True

        Returns
        -------
            A tuple of two elements:
            - The HHsearch output as a string
            - HHsearch output as a Pandas DataFrame
              See rotifer.io.hhsuite

        See also
        --------
            rotifer environment configuration

        Examples
        --------
        Load alignment in multi-FASTA format and compare it to Pfam

        >>> aln = sequence("myaln.aln")
        >>> (hhout,hhdf) = aln.hhsearch(databases=['pfam'])
        """
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer.io.hhsuite import read_hhr

        dbs = " ".join([ " -d " + os.path.join(database_path, x) for x in databases ])
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/seqaln')
            child = f'hhsearch -i {tmpdirname}/seqaln {dbs} -M 50 -cpu 18 -o {tmpdirname}/seqaln.hhr'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            hhtable = read_hhr(f'{tmpdirname}/seqaln.hhr')
            with open(f'{tmpdirname}/seqaln.hhr') as f:
                hhsearch_result = f.read()
            if view:
                from IPython.core.page import page 
                page(hhsearch_result)
            return (hhsearch_result, hhtable)

    def community(self):
        import tempfile
        import community
        import networkx as nx
        import numpy as np

        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/tt', remove_gaps=True)
            os.system(f'mmseqs easy-search {tmpdirname}/tt {tmpdirname}/tt {tmpdirname}/tt.m8 {tmpdirname}/tmp')
            df = pd.read_csv(f'{tmpdirname}/tt.m8',sep="\t", names='source target pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split())
            dfc = df
            dfc.source, dfc.target = np.where(dfc.source > dfc.target , [dfc.source, dfc.target], [dfc.target, dfc.source])
            dfc = dfc.sort_values('evalue').drop_duplicates(['source', 'target'])
            dfc['evalue_t'] = - np.log10(dfc.evalue)
            G = nx.from_pandas_edgelist(dfc[['source', 'target', 'evalue_t']].query('evalue_t >= 3'), edge_attr='evalue_t')
            partition = community.best_partition(G,weight='weight')
            c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
                {'index': 'id', 0: 'community'}, axis=1)
        return c

    def trim(self, max_perc_gaps=80, minimum_length=1):
        '''
        Remove alignment columns based on column statistics.

        Parameters
        ----------
        max_perc_gaps: integer or float
          Maximum relative frequency of gaps

        minimum_length : integer, default 1
          Minimum number of consecutive bad quality columns.
          Columns located in regions shorter than this threshold
          will not removed.

        Examples
        --------
        Remove all columns with more than 70% of gaps.
        >>> aln.trim(70)
        '''
        columns_to_keep = self.residue_frequencies.T.query('gap <= @max_perc_gaps').T.columns.to_list()
        result = self.copy()
        result.df['sequence'] = result.residues.loc[:, columns_to_keep].sum(axis=1)
        other = set(self.numerical.columns) - set(['type']) - set(list(range(1, self.get_alignment_length() + 1)))
        result.numerical = result.numerical[['type'] + columns_to_keep + list(other)]
        result.numerical.columns = ['type'] + list(range(1,len(columns_to_keep)+1)) + list(other)
        result._reset()
        return result

    def add_jpred(self, email=False, symbol=False):
        ''' 
        Function to add secondary structure from the Jpred server
        ∞ = Alpha helix
        > = Beta strand
        - = Turn
        '''
        import tarfile
        import tempfile
        import re
        from subprocess import Popen, PIPE, STDOUT, check_output
        import os
        import io

        if not email:
            from rotifer.db.ncbi import NcbiConfig
            email = NcbiConfig['email']

        result = self.copy()
        with tempfile.TemporaryDirectory() as tmpdirname:
            cd = os.getcwd()
            os.chdir(tmpdirname)
            self.to_file('seqaln')
            child = f'jpredapi submit mode=msa format=fasta email={email} file=seqaln name=rotifer'
            child = check_output(child, shell=True).decode()
            match = re.findall(f'jobid=(.+?)\s', child, re.DOTALL)[0]
            child = f'jpredapi status jobid={match} getResults=yes checkEvery=10'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            with tarfile.open(f'{match}/{match}.tar.gz', 'r:gz') as tar:
                ss_file = tar.getmember(f'{match}.jnet')
                query_file = tar.getmember(f'{match}.msf.query')
                jnet = pd.read_csv(io.BytesIO(tar.extractfile(ss_file).read()), sep=":", names = ['a', 'b'])
                query = tar.extractfile(query_file).read()
                query = re.findall('>(.+?)\n', query.decode(), re.DOTALL)[0]

            os.chdir(cd)
            jnet = jnet.iloc[0:2,:]
            jnet.b = jnet.b.str.replace(',', '')
            jnet.iloc[0,1] = jnet.query('a  == "jnetpred"').iloc[0,1].replace(',', '')
            if symbol:
                jnet.iloc[0,1] = jnet.iloc[0,1].replace('E', '>').replace('H', "=")

            a1 = pd.Series(list(result.df.query('id ==@query').iloc[0,1])).where(lambda x: x !='-').dropna().rename('seq').reset_index()
            a2 = pd.Series(list(jnet.loc[0, 'b'])).rename('jnet')
            a3 = pd.Series(list(jnet.loc[1, 'b'])).rename('conf')
            a4 = pd.Series(list(result.df.query('id == @query').iloc[0,1]))
            a5 = a1.join(a2).join(a3).set_index('index')
            a6 = pd.concat([a4,a5], axis = 1).fillna(' ').sum().reset_index()
            a6.columns = ['id', 'sequence']
            a6['type'] = 'structure prediction'
            a6.iloc[2:3,2] = 'Jpred'
            a6.iloc[3:4,2] = 'confidance'

            result.df = pd.concat([a6.iloc[2:4,:],result.df])

            return result

    def _to_df_style(self, consensus, annotations=False, remove_gaps=False):
        """TODO: Docstring for function.

        :consensus: The consensus threshold that should be used to color the aligment
        :output_file: output file name
        :annotation: List of annotations rows that should be keept in the  html file
        The annotation label should be the same as in the id seq object df columm
        :remove_gaps: Query sequence to use as model to remove the gaps, 
        it will add numbers of aminoacid suppressed in the sequence that contain the insertions.
        :returns: TODO

        """

        import sys
        import pandas as pd
        from rotifer.devel.beta.sequence import sequence
        import numpy as np
        aln = self.copy()
        if remove_gaps:
            gaps,gdf = aln.gaps_to_numbers(remove_gaps)
            gdf.index = aln.df.query('type == "sequence"').id

        aromatic = ['F','Y', 'W', 'H']
        alifatic = ['I','V','L']
        hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
        positive = ['H', 'K', 'R']
        negative = [ 'D', 'E']
        charged = positive + negative
        polar = charged + ['p','Q', 'N', 'S', 'T','C']
        alcohol = ['S','T']
        tiny = ['G', 'A', 'S']
        small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
        big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
        all_aa = ['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T']

        aa_groups_colors = {'a':[aromatic,  '#2C68F3'],
                            'l':[alifatic, '#2CF3EA'],
                            'h':[hydrophobic,  '#F3E42C90'],
                            '+':[positive,  '#2C68F3'],
                            '-':[negative,  '#F50EF195'],
                            'c':[charged,  '#38F50E'],
                            'p':[polar,  '#0EF5A150'],
                            'o':[alcohol,  '#AE5BF8'],
                            'u':[tiny,  '#EE9C0C'],
                            's':[small,  '#DA147750'],
                            'b':[big,  '#A28694'],
                            '.':[all_aa,  'white'],
                            'G':[all_aa,'white'],
                            'A':[all_aa,'white'],
                            'V':[all_aa,'white'],
                            'I':[all_aa,'white'],
                            'L':[all_aa,'white'],
                            'M':[all_aa,'white'],
                            'F':[all_aa,'white'],
                            'Y':[all_aa,'white'],
                            'W':[all_aa,'white'],
                            'H':[all_aa,'white'],
                            'C':[all_aa,'white'],
                            'P':[all_aa,'white'],
                            'K':[all_aa,'white'],
                            'R':[all_aa,'white'],
                            'D':[all_aa,'white'],
                            'E':[all_aa,'white'],
                            'Q':[all_aa,'white'],
                            'N':[all_aa,'white'],
                            'S':[all_aa,'white'],
                            'T':[all_aa,'white'],
                            ' ':[all_aa,'white'],
                            '  ':[all_aa,'white'],
                            '_':[all_aa,'white']}

        # Geting the residues tabele:
        aln_r = aln.residues
        con = pd.Series(list(aln.consensus(consensus)))
        aln_r = aln_r.set_index(aln.df.query('type == "sequence"').id)
        con.index +=1
        aln_r = pd.concat([aln_r, con.rename('consensus').to_frame().T], axis=0)
        if annotations:
            if isinstance(annotations, str):
                ann = pd.Series(list(aln.df.query('id ==@annotations').sequence.iloc[0]))
                ann.index +=1
                aln_r = pd.concat([ann.rename(annotations).to_frame().T,aln_r])
            else:
                for x in annotations:
                    ann = pd.Series(list(aln.df.query('id ==@x').sequence.iloc[0]))
                    ann.index +=1
                    aln_r = pd.concat([ann.rename(x).to_frame().T,aln_r])

        if remove_gaps:
            aln_r =  aln_r.drop(gaps,axis=1).join(gdf).sort_index(axis=1).fillna(0).astype(int, errors='ignore').astype(str).replace('0','  ')

        # Funtion that works!!!
        def highlight_aln(s):
            import numpy as np
            d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            return np.where(
                s == '  ',
                'color:white;background-color:white',
                np.where(
                    s == s.iloc[-1],
                    'color:white;background-color:black',
                    np.where(
                        s.isin(d[0]),
                        f'color:black;background-color:{d[1]}',
                        'color:black;background-color:white')))

        def highlight_consensus(s):
            import numpy as np
            d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            """TODO: Docstring for highlight_consensus.

            :arg1: TODO
            :returns: TODO

            """
            return np.where(
                s.isin(all_aa),
                'color:white;background-color:black',
                f'color:black;background-color:{d[1]}',
                )

        #Making slice index where the functions should be applied:
        # One function should be applien only in the consensus row
        # Other function should be appplied only in seq rows
        idx1 = pd.IndexSlice
        corte = idx1[idx1['consensus'],idx1[:]]
        #Getting the firs sequence (fs) row to map the slice:
        idx2 = pd.IndexSlice
        fs = aln.df.query('type == "sequence"').id.iloc[0]    
        corte2 = idx2[idx2[fs:],idx2[:]]

        headers = {
            'selector': 'th:not(.index_name)',
            'props': '''font-size: 12px;
            text-align: left;
            font-family:"Lucida Console", Monaco, monospace;
            color:black;
            background-color:white'''
        }

        if sys.version_info.minor > 8:
            df_style = aln_r.style.set_properties(**{
                'font-size': '12px',
                'font-family':'"Lucida Console", Monaco,monospace',
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide(axis='columns').apply(
                highlight_consensus, subset=corte
            ).set_table_styles(
                [headers]
            )
        else:
            df_style = aln_r.style.set_properties(**{
                'font-size': '12px',
                'font-family':'"Lucida Console", Monaco,monospace',
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide_columns().apply(
                highlight_consensus, subset=corte
            ).set_table_styles(
                [headers]
            )
        #if whant to send to latex, replace set_stick... to:to_latex(environment='longtable', convert_css=True)
        return df_style

    def to_html(self, consensus,output_file, annotations=False, remove_gaps=False, fixed_index=True):
        """TODO: Docstring for function.

        :consensus: The consensus threshold that should be used to color the aligment
        :output_file: output file name
        :annotation: List of annotations rows that should be keept in the  html file
        The annotation label should be the same as in the id seq object df columm
        :remove_gaps: Query sequence to use as model to remove the gaps, 
        it will add numbers of aminoacid suppressed in the sequence that contain the insertions.
        :fixed_index, It makes the index fixed on the html page 
        :returns: TODO

        """
        html = self._to_df_style(
            consensus,
            annotations=annotations,
            remove_gaps=remove_gaps
        )
        if fixed_index:
            html = html.set_sticky(axis="index")
            
        with open(output_file, 'w') as f:
            f.write(html.render(table_attributes='cellspacing=0, cellpadding=0'))
            return(f'{output_file} saved on the working path')


    def gaps_to_numbers(self, smodel):
        """: This function will retunr a tuple. The first element of the tuple is 
        a list containing the index of columns that contains gap on the model sequence.
        It should be used to drop the columns from the residues df.
        >>> aln.residues.drop(firs element of the tuple output, axis =1).
        
        The seccond element is a df that should be joined with the residues df after the drop:
        >>> aln.residues.join(second elemnet of tuple).sort_index(axis=1)
        :smodel: Sequence used as model to remove the gaps.

        """
        midx = self.df.query('id == @smodel').index
        gaps = self.residues.loc[midx].iloc[0].where(lambda x: x=='-').dropna()
        grouper_map = gaps.reset_index().drop(midx, axis=1).rename(
            {'index':'region'},
            axis=1
        ).eval(
            'inter = ~region.diff().fillna(1).le(1)'
        ).eval(
            'grouper = inter.cumsum()'
        )
        gap_df = grouper_map.groupby('grouper').agg(
            first = ('region','min'),
            last = ('region', 'max'),
            rsize = ('region', 'nunique'),
            rlist = ('region', list)
        )
        num_gaps = grouper_map.set_index('region')[['grouper']].join(
            (self.residues[gaps.index] != '-').T
        ).groupby('grouper').sum().T
        num_gaps.columns = gap_df['first'].to_list()
        return (list(gaps.index),num_gaps)

    def _to_df_styleTEX(self, consensus, annotations=False, remove_gaps=False):
        """TODO: Docstring for function.

        :consensus: The consensus threshold that should be used to color the aligment
        :output_file: output file name
        :annotation: List of annotations rows that should be keept in the  html file
        The annotation label should be the same as in the id seq object df columm
        :remove_gaps: Query sequence to use as model to remove the gaps, 
        it will add numbers of aminoacid suppressed in the sequence that contain the insertions.
        :returns: TODO

        """
        import sys
        import pandas as pd
        from rotifer.devel.beta.sequence import sequence
        import numpy as np
        aln = self.copy()
        if remove_gaps:
            gaps,gdf = aln.gaps_to_numbers(remove_gaps)
            gdf.index = aln.df.query('type == "sequence"').id

        aromatic = ['F','Y', 'W', 'H']
        alifatic = ['I','V','L']
        hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
        positive = ['H', 'K', 'R']
        negative = [ 'D', 'E']
        charged = positive + negative
        polar = charged + ['p','Q', 'N', 'S', 'T','C']
        alcohol = ['S','T']
        tiny = ['G', 'A', 'S']
        small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
        big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
        all_aa = ['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T']

        aa_groups_colors = {'a':[aromatic,  '#2C68F3'],
                            'l':[alifatic, '#2CF3EA'],
                            'h':[hydrophobic,  '#F3E42C90'],
                            '+':[positive,  '#2C68F3'],
                            '-':[negative,  '#F50EF195'],
                            'c':[charged,  '#38F50E'],
                            'p':[polar,  '#0EF5A150'],
                            'o':[alcohol,  '#AE5BF8'],
                            'u':[tiny,  '#EE9C0C'],
                            's':[small,  '#DA147750'],
                            'b':[big,  '#A28694'],
                            '.':[all_aa,  'white'],
                            'G':[all_aa,'white'],
                            'A':[all_aa,'white'],
                            'V':[all_aa,'white'],
                            'I':[all_aa,'white'],
                            'L':[all_aa,'white'],
                            'M':[all_aa,'white'],
                            'F':[all_aa,'white'],
                            'Y':[all_aa,'white'],
                            'W':[all_aa,'white'],
                            'H':[all_aa,'white'],
                            'C':[all_aa,'white'],
                            'P':[all_aa,'white'],
                            'K':[all_aa,'white'],
                            'R':[all_aa,'white'],
                            'D':[all_aa,'white'],
                            'E':[all_aa,'white'],
                            'Q':[all_aa,'white'],
                            'N':[all_aa,'white'],
                            'S':[all_aa,'white'],
                            'T':[all_aa,'white'],
                            ' ':[all_aa,'white'],
                            '  ':[all_aa,'white'],
                            '_':[all_aa,'white']}

        # Geting the residues tabele:
        aln_r = aln.residues
        con = pd.Series(list(aln.consensus(consensus)))
        aln_r = aln_r.set_index(aln.df.query('type == "sequence"').id)
        con.index +=1
        aln_r = pd.concat([aln_r, con.rename('consensus').to_frame().T], axis=0)
        if annotations:
            if isinstance(annotations, str):
                ann = pd.Series(list(aln.df.query('id ==@annotations').sequence.iloc[0]))
                ann.index +=1
                aln_r = pd.concat([ann.rename(annotations).to_frame().T,aln_r])
            else:
                for x in annotations:
                    ann = pd.Series(list(aln.df.query('id ==@x').sequence.iloc[0]))
                    ann.index +=1
                    aln_r = pd.concat([ann.rename(x).to_frame().T,aln_r])

        if remove_gaps:
            aln_r =  aln_r.drop(gaps,axis=1).join(gdf).sort_index(axis=1).fillna(0).astype(int, errors='ignore').astype(str).replace('0','  ')

        # Funtion that works!!!
        def highlight_aln(s):
            import numpy as np
            d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            return np.where(
                s == '  ',
                'color:white;background-color:white',
                np.where(
                    s == s.iloc[-1],
                    'color:white;background-color:black',
                    np.where(
                        s.isin(d[0]),
                        f'color:black;background-color:{d[1]}',
                        'color:black;background-color:white')))

        def highlight_consensus(s):
            import numpy as np
            d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            """TODO: Docstring for highlight_consensus.

            :arg1: TODO
            :returns: TODO

            """
            return np.where(
                s.isin(all_aa),
                'color:white;background-color:black',
                f'color:black;background-color:{d[1]}',
                )

        #Making slice index where the functions should be applied:
        # One function should be applien only in the consensus row
        # Other function should be appplied only in seq rows
        idx1 = pd.IndexSlice
        corte = idx1[idx1['consensus'],idx1[:]]
        #Getting the firs sequence (fs) row to map the slice:
        idx2 = pd.IndexSlice
        fs = aln.df.query('type == "sequence"').id.iloc[0]    
        corte2 = idx2[idx2[fs:],idx2[:]]

        if sys.version_info.minor > 8:
            df_style = aln_r.style.set_properties(**{
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide(axis='columns').apply(
                highlight_consensus, subset=corte
            )
        else:
            df_style = aln_r.style.set_properties(**{
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide_columns().apply(
                highlight_consensus, subset=corte
            )
        #if whant to send to latex, replace set_stick... to:to_latex(environment='longtable', convert_css=True)
        return df_style


    def edit(self, consensus=True, scale=True):
        """
        Search the alignment against a HMM databases using hhsearch.

        Parameters
        ----------
        databases : list of strings, default ['pfam','pdb70']
            List of HMM databases to include in the search
        database_path : string, default is ROTIFER_DATA/hhsuite
            Path to the directory where the HMM databases are stored
        view : bool, default True

        Returns
        -------
            A tuple of two elements:
            - The HHsearch output as a string
            - HHsearch output as a Pandas DataFrame
              See rotifer.io.hhsuite

        See also
        --------
            rotifer environment configuration

        Examples
        --------
        Load alignment in multi-FASTA format and compare it to Pfam

        >>> aln = sequence("myaln.aln")
        >>> (hhout,hhdf) = aln.hhsearch(databases=['pfam'])
        """
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer import GlobalConfig as gc
        import os
        os.environ["VIMMOVE"] = f'{gc["base"]}/etc/vim/vim-move'
        aln = self.copy()
        result = self.copy()
        if consensus:
            aln = aln.add_consensus(separator='=')


        alndf = aln.df.fillna('X')
        cols = list(alndf.dtypes.where(lambda x: x=='object').dropna().index) 
        if scale:
            scale = aln._scale_bar(self.get_alignment_length(), interval=10)
            alndf = pd.concat([ scale, alndf ]).fillna('X') 

        for x in cols:
            alndf[x] = alndf[x].str.pad(alndf[x].str.len().max(), side="right")
        with tempfile.TemporaryDirectory() as tmpdirname:
            alndf.to_csv(f'{tmpdirname}/seqdf.fa',index=False, sep="\t")
            print ("Open VMI to edit your file.")
            if os.system(f'nvim -u {gc["base"]}/etc/vim/init.vim {tmpdirname}/seqdf.fa') != 0:
                        raise TryNext()
            tmpdf = pd.read_csv(f'{tmpdirname}/seqdf.fa', sep="\t")
            tmpdf = tmpdf[tmpdf['type'].str.contains("sequence")]
            tmpdf = tmpdf.dropna(subset='sequence')
            result = result.filter(keep=tmpdf.id.to_list())
            result.df.sequence = result.df.id.map(tmpdf.set_index('id').sequence.dropna().to_dict())
        
        return result 





    ## Class methods

    @classmethod
    def _scale_bar(self, length, interval=10, start=1):
        # Numbers for the scale bar
        if length > interval:
            position = list(range(int(((start / interval)+1))*interval,length,interval))
            if not length % interval:
                position.append(length)
        else:
            position = [length]

        # Tick marks (vertical bars)
        scale_bar = [ f'{"|":{interval}}' for x in position ]

        # Format numbers
        scale_number = [ f'{str(x):{interval}}' for x in position ]
        if len(str(start)) < position[0] - start - 1:
            scale_number.insert(0, str(start) + " " * (position[0] - start - len(str(start))))
            scale_bar.insert(0, "|" + " " * (position[0] - start - 1))
        else:
            scale_number.insert(0, " " * (position[0] - start))
            scale_bar.insert(0, " " * (position[0] - start))

        # Prepare Series
        scale_bar = "".join(scale_bar)
        scale_number = "".join(scale_number)
        scale_number = scale_number.rstrip() + " " * (length - start - len(scale_number.rstrip()) + 1)
        scale_bar = scale_bar.rstrip() + " " * (len(scale_number) - len(scale_bar.rstrip()))
        scale_dot = "." * len(scale_number)
        scale = pd.Series([scale_number,scale_bar,scale_dot], index=['position', 'bar', 'dot'])
        scale = {'id':scale.index.tolist(), 'sequence':scale, 'length':scale.str.len(), 'type':'view'}
        scale = pd.DataFrame(scale).astype(str)
        return scale

    @classmethod
    def from_seqrecords(cls, input_data):
        """
        Build a MSA object from a list of BioPython objects.

        Parameters
        ----------
        input_data  : list of Bio.SeqRecord.SeqRecord
        """
        return cls(input_data, input_format=type(input_data[0]))

    @classmethod
    def from_string(cls, input_data, input_format='fasta'):
        '''
        This function uses BioPython to parse MSAs from strings.

        Parameters
        ----------
        input_file   : file path or open file handle
        input_format : Biopython supported file format.
                       See Bio.SeqIO and/or Bio.AlignIO.
        '''
        return cls.from_file(StringIO(input_data), input_format=input_format) 

    @classmethod
    def from_file(cls, input_file, input_format='fasta'):
        '''
        Parse multiple sequence alignment files using BioPython.

        Parameters
        ----------
        input_file   : file path or open file handle
        input_format : Biopython supported file format.
                       See Bio.SeqIO and/or Bio.AlignIO.
        '''
        return cls(input_file, input_format=input_format)

__doc__ = """
========================
Rotifer Sequence Objects
========================

This module implements classes and methods for dealing with
multiple sequence alignments, also known as MSAs.

Examples
--------
Simplest usage: visualize a (multi)FASTA file.

>>> from rotifer.devel.beta.sequence import sequence
>>> aln = sequence("alignment.aln")
>>> aln.view()

Load alignment from another file format, such as
StockHolm and show the consensus above the alignment:

>>> from rotifer.devel.beta.sequence import sequence
>>> sto = sequence("alignment.sto", format="stockholm")
>>> sto.add_consensus().view()

Parsing a FASTA alignment from an explicit string:

>>> from rotifer.devel.beta.sequence import sequence
>>> aln = sequence('>Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n')
>>> aln.add_consensus().view()
"""
