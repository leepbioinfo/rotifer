import os

# Private data
_reserved_columns = ['id','sequence','length','type']

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
      directly, i.e. using Pandas's ``append`` or ``insert`` methods
      but the user must be careful to set the value or the ``type``
      column to something other than ``sequence``, which is researved
      for the aligned sequences.

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
    freq_table : bool, default True
        Build/Update the table of per-column aminoacid frequencies.

    See also
    --------
    Bio.SeqIO   : BioPython parser for aligned/unaligned sequences
    Bio.AlignIO : BioPython parser for aligned sequences
    pandas      : library used internally to store sequences

    Examples
    --------
    Simplest usage: visualize a (multi)FASTA file.

    >>> from rotifer.devel.alpha.sequence import sequence
    >>> aln = sequence("alignment.aln")
    >>> aln.view()

    Load alignment from another file format, such as
    StockHolm and show the consensus above the alignment:

    >>> from rotifer.devel.alpha.sequence import sequence
    >>> sto = sequence("alignment.sto", format="stockholm")
    >>> sto.add_consensus().view()

    Parsing a FASTA alignment from an explicit string:

    >>> from rotifer.devel.alpha.sequence import sequence
    >>> aln = sequence(">Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n")
    >>> aln.add_consensus().view()
    """
    from IPython.core.page import page 
    def __init__(self, input_data=None, input_format='fasta', freq_table=True):
        from io import IOBase
        from io import StringIO
        from Bio import SeqIO
        self.input_format = input_format

        # Generate empty object
        if input_data is None:
            self.df = pd.DataFrame({}, columns=_reserved_columns)
            self.numerical = pd.DataFrame({}, columns=['id','type'])
            self.file_path = None
            return

        # Initialize parser for files, file handles or strings
        if isinstance(input_data,str):
            if os.path.exists(input_data):
                self.file_path = input_data
            else:
                self.file_path = 'StringIO'
                input_data = StringIO(input_data)
            input_data = SeqIO.parse(input_data, input_format)
        elif isinstance(input_data, IOBase):
            input_data = SeqIO.parse(input_data, input_format)

        # Parse each sequence and calculate frequencies
        self.df = self.__seqrecords_to_dataframe(input_data)
        if freq_table:
            self.freq_table = self._aln_freq_df(by_type=True)

    @staticmethod
    def __seqrecords_to_dataframe(data):
        import pandas as pd
        df = pd.DataFrame([ (s.id,str(s.seq),len(s.seq.replace("-","")),'sequence') for s in data ], columns=_reserved_columns)
        return df

    def __len__(self):
        return len(self.df)

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
        return self.df.sequence.str.len().max()

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
        from copy import deepcopy
        import pandas as pd
        result = deepcopy(self)

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
            result.freq_table = result._aln_freq_df(by_type=True)
        return result

    def slice(self, position):
        '''
        Method to select and concatenate one or more sets of columns.

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
        from copy import deepcopy
        import pandas as pd
        result = deepcopy(self)
        if isinstance(position[0],int):
            position = [position]

        stack = []
        for pos in position:
            pos = [*pos]
            if len(pos) == 3:
                refseq = pd.Series(list(result.df.query(f'''id == "{pos[2]}"''').sequence.values[0]))
                refseq = refseq.where(lambda x: x != '-').dropna().reset_index().rename({'index':'mapped_position'}, axis=1)
                pos[0:2] = (refseq.loc[pos[0]-1:pos[1]].mapped_position.agg(['min','max'])).tolist()
                pos[0] += 1
            stack.append(result.df.sequence.str.slice(pos[0]-1, pos[1]))
        result.df['sequence'] = pd.concat(stack, axis=1).sum(axis=1)
        result.df['length'] = result.df.sequence.str.replace('-', '').str.len()
        result.freq_table = result._aln_freq_df(by_type=True)
        return result

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
        import pandas as pd
        from ascii_graph import Pyasciigraph
        print(f'Total proteins: {len(self.df)}')
        a = self.df['length'].value_counts().to_frame().reset_index()
        a = a.sort_values('index')
        a['raw_bin'] = pd.cut(a['index'],bins,precision=0)
        a['bin'] = a.raw_bin.apply(lambda x : '{} - {}'.format(int(x.left),int(x.right)))
        test = a.groupby('bin').agg({'length':'sum'}).reset_index().apply(tuple, axis=1)
        graph = Pyasciigraph()
        for line in  graph.graph('count \t sequence size', test):
            print(line)

    def sort(self, by=['length'], ascending=True, id_list=None, tree_file=None, tree_format='newick'):
        """
        Sort alignment rows.

        Parameters
        ----------
        by : list of strings, default ['length']
            List of criteria for sorting.

            The following criteria are always supported:
            - id     : sort by sequence identifiers
            - length : sequence length
            - list   : sort by user-provided list
                See parameter ``id_list``.
            - name   : sort by sequence identifiers
            - tree   : sort by phylogeny
                See parameter ``tree_file``.
                Leaf names must match sequence identifiers.

            Additionally, user-supplied sequence annotations may
            also be used to sort alignment rows. See the section
            ``Alignment annotations``.

        ascending : bool or list of bools, default True
            If True, sort in ascending order, otherwise descending.

        id_list : list of sequence identifiers or file name
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
        Reordered copy of the original MSA

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

        """
        from copy import deepcopy
        import pandas as pd
        result = deepcopy(self)

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

            elif item in result.df.columns:
                fields.append(result.df[item])

            else:
                print(f'sequence.sort: Unsupported criteria or missing annotation ({item})', file=sys.stderr)

        # Concatenate sorting fields
        fields = pd.concat(fields, axis=1).sort_values(by=by, ascending=ascending)
        result.df = result.df.reindex(fields.index)
        return result

    def add_hhpred(self, hhpred_file = 'hhpred_file', hhpred_id = ''):
        '''
        Parse HH-suite output file and add the selected result's
        sequence and secondary structure to the MSA.
        '''
        import re
        from copy import copy, deepcopy
        import pandas as pd

        # Parser
        with open(hhpred_file, 'r') as f:
            texto = f.read()
        match = re.findall(f'No {hhpred_id}(.+?)No ', texto, re.DOTALL)[0].strip()
        target = re.findall('T Consensus.*?\nT\s(.*?)\s',match, re.MULTILINE)[0]
        query = re.findall('Q ss_pred.*?\nQ\s(.*?)\s',match, re.MULTILINE)[0]
        T_con = ''.join(re.findall('T Consensus.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        sequence_target = ''.join(re.findall(f'T\s+{target}.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        sequence_query = ''.join(re.findall(f'Q\s+{query}.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        T_ss_pred = ''.join(re.findall('T ss_pred\s+(.*?)$',match,  re.MULTILINE))
        T_ss_dssp = ''.join(re.findall('T ss_dssp\s+(.*?)$',match,  re.MULTILINE))
        Q_ss_pred = ''.join(re.findall('Q ss_pred\s+(.*?)$',match,  re.MULTILINE))
        Q_con = ''.join(re.findall('Q Consensus.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        structure_df = pd.DataFrame({'query':list(sequence_query), 'query_pred': list(Q_ss_pred), 'target':list(sequence_target), 'ss':list(T_ss_dssp)})

        # join aln to hhpred results
        result = deepcopy(self)
        aln = pd.Series(list(result.df.query('id == @query').sequence.iloc[0])).where(lambda x: x!='-').dropna().reset_index().rename({'index':'position', 0:'sequence'}, axis=1)
        s_aln = ''.join(aln.sequence).find(''.join(structure_df['query'].where(lambda x : x !='-').dropna()))
        e_aln = s_aln + len(''.join(structure_df['query'].where(lambda x : x !='-').dropna())) -1
        aln.loc[s_aln : e_aln , 'dssp'] = structure_df[structure_df['query'] != '-'].ss.to_list()
        aln.loc[s_aln : e_aln , 'sspred'] = structure_df[structure_df['query'] != '-'].query_pred.to_list()
        aln.loc[s_aln : e_aln , target] = structure_df[structure_df['query'] != '-'].target.to_list()
        aln = pd.Series(list(result.df.iloc[0].sequence)).to_frame().join(aln.set_index('position', drop=True)).fillna('-')
        aln = aln.T.apply(lambda x: ''.join(list(x)), axis=1).reset_index().rename({'index': 'id', 0: 'sequence'}, axis=1) 
        result.df =  pd.concat([aln,result.df]).iloc[2:].reset_index(drop=True)

        return result

    def add_pdb(self, pdb_id, chain_id='A', pdb_file=None):
        """
        Add structure annotation from PDB.
        """
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        from Bio import pairwise2
        import Bio.PDB.PDBList as PDBList
        from pathlib import Path
        import pandas as pd
        from copy import deepcopy
        result = deepcopy(self)
        pdb_id = pdb_id.lower()

        # Download, if needed
        if not pdb_file:
            pdb = PDBList(verbose=False).retrieve_pdb_file(pdb_code=pdb_id, file_format='pdb', pdir='./')
            p = Path(pdb)
            p.rename(p.with_suffix('.pdb'))
            pdb = pdb.replace('ent', 'pdb')
        else:
            pdb = pdb_file

        # Parse PDB file
        dssp_to_abc = {"I":"C","S":"C","H":"H","E":"E","G":"C","B":"E","T":"C","C":"C"}
        dssp = PDBParser().get_structure(pdb_id, pdb)
        dssp = DSSP(dssp[0], pdb)
        select_chain = [x for x in list(dssp.keys()) if x[0] == chain_id]
        l =  [dssp[select_chain[x]] for x in range(len(select_chain))]
        df1 = pd.DataFrame(data={'position':[x[0] for x in l],'aa': [x[1] for x in l],'structure':  [x[2] for x in l]})
        df1.structure = df1.structure.map(dssp_to_abc).fillna('-')

        # Search for pdb sequence in the aligment
        pdb_index = result.df[result.df.id.str.contains(pdb_id, case=False)].index[0]
        pdb_in_aln = result.df.loc[pdb_index].sequence.replace('-', '')
        pdb_from_pdb = ''.join(df1.aa.to_list())
        ali = pairwise2.align.localxx(pdb_in_aln, pdb_from_pdb)[0]
        ss_df = pd.Series(list(ali.seqB)).where(lambda x : x != '-').dropna().to_frame()
        ss_df['structure'] = df1.structure.to_list()
        pdb_df = pd.Series(list(ali.seqA)).rename('seq').to_frame().join(ss_df['structure']).fillna('-').query(' seq != "-"')
        to = pd.Series(list(result.df.loc[pdb_index].sequence)).where(lambda x: x !='-').dropna().rename('ung').to_frame()
        to['structure'] = pdb_df.structure.to_list()
        pdbn = f'ss_from:{pdb_id}_{chain_id}'
        pdbss = ''.join(pd.Series(list(result.df.loc[pdb_index].sequence)).to_frame().join(to).fillna('-').structure.to_list())
        result.df =  pd.concat([pd.DataFrame([[pdbn,pdbss]], columns=['id', 'sequence']),self.df])
        return result

    def to_color(self, color='fg', scale=True):
        """
        Fetch aligment as text, colored by residue class.
        """
        import pandas as pd
        df = self.df.copy()
        alignment_length = self.get_alignment_length()

        def color_res(s, cs):
            if s in 'A I L M F W V'.split():
                return cs(s, 33)
            elif s in 'K R'.split():
                return cs(s, 124)
            elif s in 'E D'.split():
                return cs(s, 127)
            elif s in 'N Q S T'.split():
                return cs(s, 34)
            elif s  == 'C':
                return cs(s, 168)
            elif s == 'G':
                return cs(s, 166)
            elif s == 'P':
                return cs(s, 178)
            elif s in 'H Y'.split():
                return cs(s, 37)
            else:
                 return s

        def color_bg(s, color = ''):
            '''
            s: String
            '''
            if color:
                color = f'48;5;{color}'
            return f'\033[{color}m{s}\033[m'

        def color_fg(s, color = ''):
            if color:
                color = f'38;5;{color}'
            return f'\033[{color}m{s}\033[m'

        color_switch = {'background':color_bg, 'bg':color_bg, 'foreground':color_fg, 'fg':color_fg}
        df['colored'] = df['sequence'].map(lambda x: ''.join([color_res(y, color_switch[color]) for y in x]))
        if scale:
            scale_number = list(range(0,alignment_length,10))
            scale_number[0] = 1
            scale_number = "".join([ f'{str(x):10}' for x in scale_number ])
            scale_number = scale_number.rstrip() + " " * (alignment_length - len(scale_number.rstrip()))
            scale_bar = "".join([ f'{"|":9}' for x in range(0,alignment_length,10) ])
            scale_dot = "".join([ "." for x in range(0,alignment_length) ])
            scaled = pd.concat([pd.Series([scale_number,scale_bar,scale_dot], index=['position', 'bar', 'dot']), df.set_index('id').colored]) 
            return scaled.str.ljust(scaled.str.len().max())
        else:
            return df.colored.str.ljust(df.colored.str.len().max())

    def _aln_freq_df(self, by_type=False):
        from copy import copy, deepcopy
        import pandas as pd
        result = deepcopy(self)
        freq_df = result.df.sequence.str.split('', expand=True).iloc[:,1:-1].apply(pd.value_counts).fillna(0).astype(int)/len(result.df)*100 
        freq_df.rename({'-':'gap','.':'gap','?':'X'}, inplace=True)
        if not by_type:
            return  freq_df

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

    def consensus(self, cons):
        from copy import copy, deepcopy
        import pandas as pd
        '''
        Generate consensus lines for the alignment object
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
        result = deepcopy(self.freq_table)
        result.rename({'gap':'.'}, inplace=True)
        result = pd.concat([result, pd.DataFrame(columns=result.columns, index=['.']).fillna(101)])
        freq = result.melt(ignore_index=False).reset_index().rename({'index':'aa', 'variable':'position', 'value':'freq'}, axis=1)
        freq['ranking'] = freq.aa.map(aa_type_dict)
        freq = freq.sort_values(['position', 'ranking'], na_position='first').query(f'freq >={cons}').drop_duplicates(subset='position')

        return ''.join(freq.aa.to_list())

    def add_consensus(self, consensus=(50, 60, 70, 80, 90)):
        from copy import deepcopy
        import pandas as pd
        result = deepcopy(self)
        for x in consensus:
            cx = self.consensus(x)
            result.df = pd.concat([pd.DataFrame([[f'r_seq consensus{x}%',cx]], columns=['id', 'sequence']),result.df])
        return result

    def to_seqrecords(self, annotations=None, remove_gaps=False):
        """
        Convert the MSA to a list of BioPython's SeqRecord objects.
        """
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # Access and filter internal dataframe
        df = self.df
        if not annotations:
            df = df[~df.id.str.startswith('r_seq')]

        # Convert each row to a SeqRecord
        result = []
        for row in df.values:
            s = row[1]
            if remove_gaps:
                s = s.replace("-","")
            result.append(SeqRecord(id=row[0], seq=Seq(s)))
        return result

    def to_file(self, file_path=None, output_format='fasta', annotations=None, remove_gaps=False):
        from Bio import SeqIO
        return SeqIO.write(self.to_seqrecords(annotations=annotations, remove_gaps=remove_gaps), file_path, output_format)

    def to_string(self, output_format='fasta', annotations=None, remove_gaps=False):
        from Bio import SeqIO
        from io import StringIO
        sio = StringIO("")
        SeqIO.write(self.to_seqrecords(annotations=annotations, remove_gaps=remove_gaps), sio, output_format)
        return sio.getvalue()

    def view(self, color=True, annotations=False):
        from IPython.core.page import page
        if color:
            page(self.to_color().__repr__())
        else:
            page(self.df.set_index('id').sequence.__repr__())

    def realign(self,fast=False, cpu=10):
        from subprocess import Popen, PIPE, STDOUT
        seq_string = self.to_string().encode()
        if fast:
            child = Popen(f'cat|mafft  --thread {cpu} -' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        else:
            child = Popen(f'cat|mafft  --maxiterate 1000 --localpair --thread {cpu} -' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        result = self.from_string(child[0].decode("utf-8"), input_format = 'fasta')
        result.file_path = 'from realign function'
        #result.freq_table = result._aln_freq_df(by_type=True) 
        return result

    def hhsearch(self, databases=['pdb70','pfam'], database_path=os.path.join(os.environ['ROTIFER_DATA'],"hhsuite")):
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer.io.hhsuite import parse_hhr

        dbs = " ".join([ " -d " + os.path.join(database_path, x) for x in databases ])
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/seqaln')
            child = f'hhsearch -i {tmpdirname}/seqaln {dbs} -M 50 -cpu 18 -o {tmpdirname}/seqaln.hhr'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            hhtable = parse_hhr(f'{tmpdirname}/seqaln')
            with open(f'{tmpdirname}/seqaln.hhr') as f:
                hhsearch_result = f.read()
            return (hhsearch_result, hhtable)

    def community(self):
        import tempfile
        import os
        import pandas as pd
        import community
        import networkx as nx
        import numpy as np

        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/tt')
            os.system(f'mmseqs easy-search {tmpdirname}/tt {tmpdirname}/tt {tmpdirname}/tt.m8 {tmpdirname}/tmp')
            df = pd.read_csv(f'{tmpdirname}/tt.m8',sep="\t", names='source target pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split())
            dfc = df
            dfc.source, dfc.target = np.where(dfc.source > dfc.target , [dfc.source, dfc.target], [dfc.target, dfc.source])
            dfc = dfc.sort_values('evalue').drop_duplicates(['source', 'target'])
            dfc['evalue_t'] = - np.log10(dfc.evalue)
            G = nx.from_pandas_edgelist(dfc[['source', 'target', 'evalue_t']].query('evalue_t >= 3'), edge_attr='evalue_t')
            partition = community.best_partition(G,weight='weight')
            c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
                {'index': 'c80e3', 0: 'community'}, axis=1)
        return c

    def add_seq(self, seq_to_add, cpu=12, fast=False):
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        import pandas as pd

        # Dump input sequences
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Save current sequences to temporary file
            self.to_file(f'{tmpdirname}/seqaln') 

            # Save new sequences to temporary file
            from Bio.SeqRecord import SeqRecord
            if isinstance(seq_to_add[0],SeqRecord):
                from Bio import SeqIO
                SeqIO.write(seq_to_add, f'{tmpdirname}/acc.fa', "fasta")
            else:
                pd.Series(seq_to_add).to_csv(f'{tmpdirname}/acc', index=None, header=None)
                Popen(f'pfetch {tmpdirname}/acc > {tmpdirname}/acc.fa' , stdout=PIPE,shell=True).communicate() 

            # Run MAFFT
            if fast:
                child = Popen(f'mafft  --thread {cpu} --add {tmpdirname}/acc.fa {tmpdirname}/seqaln', stdout=PIPE,shell=True).communicate()
            else:
                child = Popen(f'mafft  --maxiterate 1000 --localpair --thread {cpu} --add {tmpdirname}/acc.fa  {tmpdirname}/seqaln', stdout=PIPE,shell=True).communicate()

        # Load output alignment
        result = self.from_string(child[0].decode("utf-8"), input_format = 'fasta')
        result.file_path = 'from add_seq function'
        #result.freq_table = result._aln_freq_df(by_type=True) 
        return result

    def add_cluster(self, coverage=0.8, identity=0.7):
        '''
        Add or update MMseqs2 clustering data for all sequences.

        Parameters
        ----------
        coverage : float
            Minimum coverage required for both sequences in each
            pairwise alignment
        identity : float
            Minimum percentage identity required for both sequences
            in each pairwise alignment

        Returns
        -------
        A new MSA object
        '''
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from copy import deepcopy
        import pandas as pd
        import os

        result = deepcopy(self)
        path = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            result.to_file(f'{tmpdirname}/seqaln', remove_gaps=True)
            Popen(f'mmseqs easy-cluster {tmpdirname}/seqaln nr --min-seq-id {identity} -c {coverage} tmp', stdout=PIPE,shell=True).communicate()
            d = pd.read_csv(f'{tmpdirname}/nr_cluster.tsv', sep="\t", names=['cluster', 'pid']).set_index('pid').cluster.to_dict()
            result.df[f'c{int(float(coverage)*100)}i{int(float(identity)*100)}'] = result.df.id.replace(d)
            os.chdir(path)

        return result

    ## Class methods

    @classmethod
    def from_seqrecords(cls, input_data, freq_table=True):
        """
        Build a MSA object from a list of BioPython objects.

        Parameters
        ----------
            input_data : list of Bio.SeqRecord.SeqRecord
            freq_table : bool, default True
                         Calculate amino acid frequency table
        """
        return cls(input_data, input_format=type(input_data[0]), freq_table=freq_table)

    @classmethod
    def from_string(cls, input_data, input_format='fasta', freq_table=True):
        '''
        This function uses BioPython to parse MSAs from strings.

        Parameters
        ----------
            input_file   : file path or open file handle
            input_format : Biopython supported file format.
                           See Bio.SeqIO and/or Bio.AlignIO.
            freq_table   : bool, default True
                           Calculate amino acid frequency table
        '''
        from io import StringIO
        return cls.from_file(StringIO(input_data), input_format=input_format, freq_table=freq_table) 

    @classmethod
    def from_file(cls, input_file, input_format='fasta', freq_table=True):
        '''
        Parse multiple sequence alignment files using BioPython.

        Parameters
        ----------
            input_file   : file path or open file handle
            input_format : Biopython supported file format.
                           See Bio.SeqIO and/or Bio.AlignIO.
            freq_table   : bool, default True
                           Calculate amino acid frequency table
        '''
        return cls(input_file, input_format=input_format, freq_table=freq_table)

__doc__ = """
========================
Rotifer Sequence Objects
========================

This module implements classes and methods for dealing with
multiple sequence alignments, also known as MSAs.

Examples
--------
Simplest usage: visualize a (multi)FASTA file.

>>> from rotifer.devel.alpha.sequence import sequence
>>> aln = sequence("alignment.aln")
>>> aln.view()

Load alignment from another file format, such as
StockHolm and show the consensus above the alignment:

>>> from rotifer.devel.alpha.sequence import sequence
>>> sto = sequence("alignment.sto", format="stockholm")
>>> sto.add_consensus().view()

Parsing a FASTA alignment from an explicit string:

>>> from rotifer.devel.alpha.sequence import sequence
>>> aln = sequence('>Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n')
>>> aln.add_consensus().view()
"""
