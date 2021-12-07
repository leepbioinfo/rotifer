class sequence:
    import pandas as pd
    from Bio import SeqIO
    from io import StringIO

    def __init__(self,file_path, seq_type = 'fasta', input_is_text=False, freq_table=True):
        import pandas as pd
        from Bio import SeqIO
        from io import StringIO
        '''
        Function that uses the SeqIO parser to transform sequences from multiples
        format  in  in a pandas DataFrame
        '''
        id_list, seq_list = [], []
        if input_is_text:
            fasta_sequences = SeqIO.parse(StringIO(file_path), seq_type)
        else:
            fasta_sequences = SeqIO.parse(file_path, seq_type)
        for fasta in fasta_sequences:
            id_list.append(fasta.id)
            seq_list.append(str(fasta.seq))
        data = {'id': id_list, 'sequence': seq_list}
        df = pd.DataFrame(data=data)
        df['len'] = df.sequence.str.replace('-', '').str.len()
        self.df = df
        self.file_path = file_path
        self.input_type = seq_type
        self._fromText = input_is_text
        if freq_table:
            self.freq_table = self._aln_freq_df(by_type=True)

    def filter_df(
        self,
        by=None,
        lower_limit=0,
        upper_limit=0,
        id_list=[],
        reg_ex=None
    ):

        '''
        Function to filter the df aligment
        '''
        from copy import copy, deepcopy
        import pandas as pd
        result = deepcopy(self)
        if by == 'id_list':
            result.df =  result.df.query('id in @id_list')
            return result
        elif by == 'size':
            result.df =  result.df.query('len >= @lower_limit and len <= @upper_limit')
            return result

        elif reg_ex:
            result.df = result.df.query('sequence.str.replace("-", "").str.contains(@reg_ex, regex=True)')
            return result

    def slice_sequence(self, position, query=False):
        from copy import copy, deepcopy
        import pandas as pd
        result = deepcopy(self)
        '''
        Slice the aligment by a given pair of cordinates a tuple (start and end coordinates)
        If a query is given it will slice the aligment based on the position the query cordinantes 
        on the aliment
        '''
        if not query:
            result.df.sequence = result.df.sequence.str.slice(position[0], position[1])
            result.df['len'] = result.df.sequence.str.replace('-', '').str.len()
            return result
        if query:
            cord = pd.Series(list(result.df[result.df.id == query].sequence.values[0])).where(
                lambda x: x != '-'
            ).dropna().reset_index().rename(
                {'index':'maped_position'}, axis=1
            ).loc[position[0]:position[1]-1].maped_position.describe()
            result.df.sequence = result.df.sequence.str.split('', expand=True).loc[:, cord['min']   :cord['max']+1].values.sum(axis=1)
            result.df['len'] = result.df.sequence.str.replace('-', '').str.len()
            return result

    def to_hist(self,column = 'len', bins=10):
        import pandas as pd
        from ascii_graph import Pyasciigraph
        print(f'Total proteins: {len(self.df)}')
        a = self.df['len'].value_counts().to_frame().reset_index()
        a = a.sort_values('index')
        a['raw_bin'] = pd.cut(a['index'],bins,precision=0)
        a['bin'] = a.raw_bin.apply(lambda x : '{} - {}'.format(int(x.left),int(x.right)))
        test = a.groupby('bin').agg({'len':'sum'}).reset_index().apply(tuple, axis=1)
        graph = Pyasciigraph()
        for line in  graph.graph('count \t sequence size', test):
            print(line)

    def order_df(self, by='size', tree_file='tree_file'):
        from copy import copy, deepcopy
        import pandas as pd
        result = deepcopy(self)

        if by == 'tree':
            tree = Tree(tree_file)
            leaves_name = [x.name for x in tree.get_leaves()]
            leaves_name = [x.strip("'") for x in leaves_name] 
            to = pd.DataFrame(leaves_name,columns=['ID'])
            result.df = to.merge(result.df)
            return  result

        elif by == 'file':
            to =pd.read_csv(order_file, names=['ID'])
            result.df =   to.merge(result.df)
            return  result

        elif by == 'size':
            result.df =  result.df.sort_values('len', ascending=False)
            return  result

        elif by == 'name':
            result.df =  result.df.sort_values('ID', ascending=True, inplace=True)
            return  result


    def add_hhpred(self, hhpred_file = 'hhpred_file', hhpred_id = ''):
        import re
        from copy import copy, deepcopy
        import pandas as pd
        '''
        Parser HHPRED file and add the result to a sequence df
        '''
        with open(hhpred_file, 'r') as f:
            texto = f.read()

        result = deepcopy(self)
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
        aln = pd.Series(list(result.df.query('id == @query').sequence.iloc[0])).where(lambda x : x!='-').dropna().reset_index().rename({'index':'position', 0:'sequence'}, axis=1)
        s_aln = ''.join(aln.sequence).find(''.join(structure_df['query'].where(lambda x : x !='-').dropna()))
        e_aln = s_aln + len(''.join(structure_df['query'].where(lambda x : x !='-').dropna())) -1
        aln.loc[s_aln : e_aln , 'dssp'] = structure_df[structure_df['query'] != '-'].ss.to_list()
        aln.loc[s_aln : e_aln , 'sspred'] = structure_df[structure_df['query'] != '-'].query_pred.to_list()
        aln.loc[s_aln : e_aln , target] = structure_df[structure_df['query'] != '-'].target.to_list()
        aln = pd.Series(list(result.df.iloc[0].sequence)).to_frame().join(aln.set_index('position', drop=True)).fillna('-')
        aln = aln.T.apply(lambda x: ''.join(list(x)), axis=1).reset_index().rename({'index': 'id', 0: 'sequence'}, axis=1) 
        result.df =  pd.concat([aln,result.df]).iloc[2:].reset_index(drop=True)

        return result

    def add_pdb_to_aln(self, pdb_name, pdb_file=None, chain_id='A'):
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        from Bio import pairwise2
        import Bio.PDB.PDBList as PDBList
        from pathlib import Path
        import pandas as pd
        from copy import copy, deepcopy

        result = deepcopy(self)
        pdb_name = pdb_name.lower()
        a = PDBList(verbose=False)
        if not pdb_file:
            pdb = a.retrieve_pdb_file(pdb_code=pdb_name, file_format='pdb', pdir='./')
            p = Path(pdb)
            p.rename(p.with_suffix('.pdb'))
            pdb = pdb.replace('ent', 'pdb')
        else:
            pdb = pdb_file

        p = PDBParser()
        structure = p.get_structure(pdb_name, pdb)
        model = structure[0]
        dssp = DSSP(model, pdb)
        dssp_to_abc = {"I" : "C",
                   "S" : "C",
                   "H" : "H",
                   "E" : "E",
                   "G" : "C",
                   "B" : "E",
                   "T" : "C",
                   "C" : "C"}
        a_keys = list(dssp.keys())
        select_chain = [x for x in a_keys if x[0] == chain_id]
        l =  [dssp[select_chain[x]] for x in range(len(select_chain))]
        df1 = pd.DataFrame(data={'position':[x[0] for x in l],'aa': [x[1] for x in l],'structure':  [x[2] for x in l]})
        df1.structure = df1.structure.map(dssp_to_abc).fillna('-')
        # Search for pdb sequence in the aligment
        pdb_index = result.df[result.df.id.str.contains(pdb_name, case=False)].index[0]
        pdb_sequence = pd.Series(list(result.df.sequence[pdb_index])).where(lambda x: x !='-').dropna().rename('ung').to_frame()
        pdb_in_aln = result.df.loc[pdb_index].sequence.replace('-', '')
        pdb_from_pdb = ''.join(df1.aa.to_list())
        ali = pairwise2.align.localxx(pdb_in_aln, pdb_from_pdb)[0]
        ss_df = pd.Series(list(ali.seqB)).where(lambda x : x != '-').dropna().to_frame()
        ss_df['structure'] = df1.structure.to_list()
        pdb_df = pd.Series(list(ali.seqA)).rename('seq').to_frame().join(ss_df['structure']).fillna('-').query(' seq != "-"')
        to = pd.Series(list(result.df.loc[pdb_index].sequence)).where(lambda x: x !='-').dropna().rename('ung').to_frame()
        to['structure'] = pdb_df.structure.to_list()
        pdbn = f'ss_from:{pdb_name}_{chain_id}'
        pdbss = ''.join(pd.Series(list(result.df.loc[pdb_index].sequence)).to_frame().join(to).fillna('-').structure.to_list())
        result.df =  pd.concat([pd.DataFrame([[pdbn,pdbss]], columns=['id', 'sequence']),self.df])
        return result

    def to_color(self, color='fg', scale=True):
        import pandas as pd
        '''
        Function to output the aligment colored by residues characteristics:
            '''
        df = self.df.copy()
        scale_size = len(df['sequence'].values[0])
        scale_dot = ''
        scale_number = ''
        scale_bar = ''
        for x in range(0, scale_size):
            if x == 0:
                scale_number = str(1)
                scale_bar += '|'
                scale_dot += '.'
            if x % 10 == 0:
                scale_number =scale_number[:-1]
                scale_number += str(x)
                scale_number += ' '*(10-len(str(x+1)))
                scale_number += ' '
                scale_bar = scale_bar [:-1]
                scale_bar += '|'
                scale_bar += ' '
                scale_dot += '.'
            else:
                scale_dot += '.'
                scale_bar += ' '

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

        color_switch = {'background':color_bg,
                  'bg':color_bg,
                  'foreground': color_fg,
                  'fg': color_fg}

        df['colored'] = df['sequence'].map(lambda x: ''.join([color_res(y, color_switch[color]) for y in x]))
        
        if scale: 
            color_scaled = pd.concat([pd.Series([scale_number,scale_bar,scale_dot], index=['position', 'bar', 'dot']),df.set_index('id').colored]) 

            return color_scaled.str.ljust(color_scaled.str.len().max())
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

        aa_type = [ aromatic, alifatic, hydrophobic, positive, negative, charged, polar, alcohol, tiny, small, big, all_aa]
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

        for x in aa_type_names.keys():
            freq_df = pd.concat([freq_df,pd.DataFrame({aa_type_names[x][1]:freq_df.loc[aa_type_names[x][0]].sum()}).T])
        #freq_df = pd.concat([freq_df,pd.DataFrame(columns=freq_df.columns, index=['.']).fillna(101)])
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
        result = deepcopy(self)
        result.rename({'gap':'.'}, inplace=True)
        result = pd.concat([result, pd.DataFrame(columns=result.columns, index=['.']).fillna(101)])
        freq = result.freq_table.melt(ignore_index=False).reset_index().rename({'index':'aa', 'variable':'position', 'value':'freq'}, axis=1)
        freq['ranking'] = freq.aa.map(aa_type_dict)
        freq = freq.sort_values(['position', 'ranking'], na_position='first').query(f'freq >={cons}').drop_duplicates(subset='position')

        return ''.join(freq.aa.to_list())

    def add_consensus(self, consensus=(50, 60, 70, 80, 90)):
        from copy import copy, deepcopy
        import pandas as pd
        result = deepcopy(self)
        for x in consensus:
            cx = self.consensus(x)
            result.df = pd.concat([pd.DataFrame([[x,cx]], columns=['id', 'sequence']),result.df])
        return result

    def to_file(self, file_path=None, out_format='fasta'):
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        return SeqIO.write([ SeqRecord(id=x[0], seq=Seq(x[1])) for x in self.df.values ], file_path, out_format)
    
    def to_string(self, out_format='fasta'):
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from io import StringIO
        sio = StringIO("")
        SeqIO.write([ SeqRecord(id=x[0], seq=Seq(x[1])) for x in self.df.values ], sio, out_format)
        return sio.getvalue()

    def view(self, color=True):
        from IPython.core.page import page
        if color:
            page(self.to_color().__repr__())
        else:
            page(self.__repr__())
