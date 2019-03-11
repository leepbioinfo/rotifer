#!/usr/bin/env python3

import os
import sys
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkstemp
import rotifer.core.functions as rcf
import rotifer.table.table as tb
import pandas as pd
from pandas.compat import StringIO
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns
import matplotlib.pyplot as plt
from rotifer.core.functions import loadClasses

from Bio import SeqIO
from io import StringIO as iStringIO

__version__ = 0.10
__authors__ = 'Gilberto Kaihami; Gianlucca Nicastro'

class MsaTable(pd.DataFrame):
    '''
    This class represents as a table (pandas DataFrame) all alignments
    It is possible to use all power of numpy/pandas, like loc and iloc
    Also added some custom functions which includes:
        - slice
        - hydrophobic
        - plot_hydrophobic
        - consensus
        - conservation (including:
                            - Shannon entropy (classic)
                            - Jensen Shannon Divergence
                            - Property Entropy
                              )

    A sequence alignment need at least 2 features:
        - sequence
        - header/id
    These features are represented as metadata. If the user give this information
    it is possible to work with an n x m dimensional array, keeping all additional information (e.g. taxonomy).

    This can be a major advantage since later it is possible to get some useful information from fragments of the df
    without compromissing the overall performance.

    For example given the following df:
    +------+----------+-------+
    | name | sequence | clade |
    | A    | PCK-ACG- | 1     |
    | B    | PCK-ACG- | 1     |
    | C    | ACT-RCG- | 2     |
    | D    | ACT-RCG- | 2     |
    +------+----------+-------+
    metadata = {'ID': 'name',
                'sequence':'sequence'}

    In the former example it is clear the advantages of using a dataframe.

    We can fast and simple calculate the conservation (entropy) for each clade:

    se_clade_1 = df[df['clade'] == 1].conservation()
    se_clade_2 = df[df['clade'] == 2].conservation()

    Since conservation method returns a numpy array it is possible to:
    delta_se = se_clade_1 - se_clade_2
    returning the difference between the entropy for each amino acid in each clade.

    TODO: Still need to add new methods for consevation, improve documentation and check if all math is right.
    TODO: The conservation algo can be improved.
    TODO: Insert metadata access when needed

    MAYBE: Add clustering?

    # Notes:
    # The hydrophobic graph and calculation ideia were derived from scikit-bio (https://github.com/biocore/scikit-bio)
    # The conservation measurements code were adapted from:
    #   Capra JA and Singh M.
    #   Predicting functionally important residues from sequence
    #   conservation. Bioinformatics. 23(15): 1875-1882, 2007

    '''
    _metadata = ['alignment',
                 'verbose',
                 'log_file',
                 'format',
                 'metadata',
                 'matrix',
                 'a',
                 'hydrophobic_matrix']

    def __init__(self, *args, **kwargs):
        self.alignment = kwargs.pop('alignment', None)
        self.verbose = kwargs.pop('verbose', None)
        self.log_file = kwargs.pop('log_file', None)
        self.format = kwargs.pop('format', None)
        self.metadata = kwargs.pop('metadata', None)
        self.matrix = ''
        self.a = ''
        self.hydrophobic_matrix = ''

        super(MsaTable, self).__init__(*args, **kwargs)

    # Pandas contructor
    @property
    def _constructor(self):
        return MsaTable

    def slice(self, seq_col = '', start = 0, end = 0, **kwargs):
        '''
        Slice an alignment
        '''
        df = self
        if 'sequence' in self.metadata.keys():
            seq_col = self.metadata['sequence']

        else:
            seq_col = seq_col

        return df[seq_col].str[start:end]

    def seq2matrix(self, seq_col = '',df = '', no_self = False,  **kwargs):
        if isinstance(df, str):
            df = self

        if 'sequence' in self.metadata.keys():
            seq_col = self.metadata['sequence']
        else:
            seq_col = seq_col
        if no_self:
            return df[seq_col].map(lambda x: list(x)).apply(pd.Series)

        else:
            self.matrix = df[seq_col].map(lambda x: list(x)).apply(pd.Series)

    def hydrophobic(self, **kwargs):
        '''
        # Data derived from AAIndex:
        # # http://www.genome.jp/dbget-bin/www_bget?aaindex:ARGP820101
        '''
        import numpy as np
        from collections import defaultdict
        self._hydrophobicity_idx = defaultdict(lambda: np.nan)
        self._hydrophobicity_idx.update({'A': 0.61, 'L': 1.53, 'R': 0.60, 'K': 1.15, 'N': 0.06, 'M': 1.18,
                                  'D': 0.46, 'F': 2.02, 'C': 1.07, 'P': 1.95, 'Q': 0., 'S': 0.05,
                                  'E': 0.47, 'T': 0.05, 'G': 0.07, 'W': 2.65, 'H': 0.61, 'Y': 1.88,
                                  'I': 2.22, 'V': 1.32, '-': 1.15 })
        self._hydrophobicity_labels=['Hydrophilic', 'Medium', 'Hydrophobic']

        if not isinstance(self.matrix, pd.DataFrame):
            self.seq2matrix(**kwargs)

        else:
            pass

        self.hydrophobic_matrix =  self.matrix.apply(lambda x: [self._hydrophobicity_idx[y] for y in x])

        self.hydrophobic_matrix = self.hydrophobic_matrix.fillna(1.15)

    def plot_hydrophobic(self, figsize = (15,10),
                         cmap = sns.diverging_palette(240, 10, as_cmap=True),
                         dpi = 300,
                         cbar_orientation = 'horizontal',
                         cbar_fraction = 0.0046,
                         cbar_pad = 0.04,
                         consensus_level = [0.9],
                         slice = ()
                         ):
        from matplotlib import colors, cm
        self.hydrophobic()

        df = self

        fig, ax = plt.subplots(figsize = figsize, dpi = dpi)

        if slice:
            matrix = self.hydrophobic_matrix.loc[:, slice[0]:slice[1]]
        else:
            slice = (0, self.hydrophobic_matrix.shape[1])
            matrix = self.hydrophobic_matrix

        cax = ax.imshow(matrix,  cmap=cmap)

        values = list(self._hydrophobicity_idx.values())

        ax.set_xticks(range(self.matrix.loc[:, slice[0]:slice[1]].shape[1]))

        # Need to fix this

        ax.set_yticklabels(df.ID.values, size=7)

        ax.set_xticklabels(self.consensus([0.9])[0.9][slice[0]:slice[1]], size=7)

        cbar = fig.colorbar(cax,
                            ticks=[min(values),
                                  np.nanmedian(values),
                                  max(values)],
                            orientation= cbar_orientation,
                            fraction=cbar_fraction, pad = cbar_pad
                            )

        if cbar_orientation == 'horizontal':
            cbar.ax.set_xticklabels(self._hydrophobicity_labels) # horizontal colorbar
        if cbar_orientation == 'vertical':
            cbar.ax.set_yticklabels(self._hydrophobicity_labels) # vertical colorbar

        fig.tight_layout()

        return ax

    # def aa_frequency2(self, gap = True, **kwargs):
    #     dc_groups = {
    #         'a' : ['F','Y', 'W', 'H'],
    #         'l' : ['I','V','L'],
    #         'h' : ['F', 'Y', 'W', 'H' 'A', 'C', 'F', 'M', 'W', 'Y'],
    #         '+' : ['H', 'K', 'R'],
    #         '-' : [ 'D', 'E'],
    #         'c' : ['H', 'K', 'R', 'D', 'E'],
    #         'p' : ['H', 'K', 'R', 'D', 'E','Q', 'N', 'S', 'T','C'],
    #         'o' : ['S','T'],
    #         'u' : ['G', 'A', 'S'],
    #         's' : ['G', 'A', 'S','V', 'T', 'D', 'N', 'P', 'C'],
    #         'b' : ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
    #     }
    #
    #     all_aa = ['G','A','V','I','L','M','F','Y',
    #               'W','H','C','P','K','R','D','E',
    #               'Q','N','S','T', '-']
    #
    #     if not isinstance(self.matrix, pd.DataFrame):
    #         self.seq2matrix(**kwargs)
    #
    #     else:
    #         pass
    #
    #     residue_matrix = pd.DataFrame(index=all_aa)
    #
    #     group_matrix=pd.DataFrame(index=['a', 'l', 'h', '+',
    #                                  '-', 'c', 'p', 'o',
    #                                  'u', 's', 'b'])
    #
    #     a = self.matrix.apply(pd.Series.value_counts)
    #     if gap:
    #         residue_matrix = residue_matrix.join(a)
    #         residue_matrix = residue_matrix.rename(index={'-': '_'})
    #
    #     else:
    #         residue_matrix = residue_matrix.join(a)
    #         residue_matrix = residue_matrix.rename(index={'-': '.'})
    #
    #     b = self.matrix.copy()
    #     b = b.values
    #
    #     for k,v in dc_groups.items():
    #         for aa in v:
    #             b[b == aa] = k
    #             b = np.where(b == aa, k, b)
    #
    #     c =  pd.DataFrame( data = b)
    #
    #     group_matrix = group_matrix.join(c.apply(pd.Series.value_counts))
    #
    #     # Improve names
    #     self.a = residue_matrix/(self.matrix.shape[0])
    #     self.c = group_matrix/(self.matrix.shape[0])
    #
    def aa_frequency(self, gap = True, **kwargs):

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

        all_aa = ['G','A','V','I','L','M','F','Y',
                  'W','H','C','P','K','R','D','E',
                  'Q','N','S','T', '-']

        self.group_size = pd.DataFrame.from_dict({'a': len(aromatic),
                                             'l': len(alifatic),
                                             'h': len(hydrophobic),
                                             '+': len(positive),
                                             '-': len(negative),
                                             'c': len(charged),
                                             'p': len(polar),
                                             'o': len(alcohol),
                                             'u': len(tiny),
                                             's': len(small),
                                             'b': len(big),
                                            }, orient='index').rename(
                                                {0:'size'}, axis=1)

        self.group_ranking = {k:v for k,v in  enumerate(['a', 'l', 'h', '+',
                                          '-', 'c', 'p', 'o',
                                          'u', 's', 'b'])}

        if not isinstance(self.matrix, pd.DataFrame):
            self.seq2matrix(**kwargs)

        else:
            pass

        residue_matrix = pd.DataFrame(index=all_aa)

        group_matrix=pd.DataFrame(index=['a', 'l', 'h', '+',
                                     '-', 'c', 'p', 'o',
                                     'u', 's', 'b'])

        for x in range(len(self.matrix.T)):
            a = self.matrix[x].value_counts().to_frame()

            if gap:
                residue_matrix = residue_matrix.join(a.rename(index={'-': '_'}))
                residue_matrix = residue_matrix.rename(index={'-': '_'})
            else:
                residue_matrix = residue_matrix.join(a.rename(index={'-': '.'}))
                residue_matrix = residue_matrix.rename(index={'-': '.'})

            d = {
                'a':np.where(self.matrix[x].isin(aromatic),1,0).sum(),
                'l':np.where(self.matrix[x].isin(alifatic),1,0).sum(),
                'h':np.where(self.matrix[x].isin(hydrophobic),1,0).sum(),
                '+':np.where(self.matrix[x].isin(positive),1,0).sum(),
                '-':np.where(self.matrix[x].isin(negative),1,0).sum(),
                'c':np.where(self.matrix[x].isin(charged),1,0).sum(),
                'p':np.where(self.matrix[x].isin(polar),1,0).sum(),
                'o':np.where(self.matrix[x].isin(alcohol),1,0).sum(),
                'u':np.where(self.matrix[x].isin(tiny),1,0).sum(),
                's':np.where(self.matrix[x].isin(small),1,0).sum(),
                'b':np.where(self.matrix[x].isin(big),1,0).sum(),
                }

            group_matrix = group_matrix.join(
                pd.DataFrame.from_dict(d,orient='index').rename({0:x}, axis=1)
                    )

        # Improve names
        self.a = residue_matrix/len(self.matrix)
        self.c = group_matrix/len(self.matrix)
        self.c = self.c.reset_index(drop = True)

    def consensus(self, thresholds,
                  output = 'dict',
                  **kwargs):
        '''
        output: dict/dataframe/matrix
        '''

        if isinstance(self.a, pd.DataFrame):
            pass
        else:
            # print('Running')
            self.aa_frequency(**kwargs)

        dc = {}
        if not isinstance(thresholds, (list, tuple)):
            thresholds = [thresholds]
        for threshold in thresholds:
            residue_matrix = self.a.fillna(0)

            residue_res=  residue_matrix[residue_matrix >= threshold].idxmax().dropna()

            group_matrix = self.c.fillna(0)

            group_res = (group_matrix[group_matrix >= threshold])
            group_res = group_res.apply(lambda x: [idx for idx, v in enumerate(x) if pd.notnull(v) ]
                            ).apply(lambda x: min(x) if x else np.nan).dropna()

            # self.debug = group_res

            group_res = group_res.map(self.group_ranking)

            e = pd.DataFrame(index = [x for x in range(0, self.matrix.shape[1])])

            e = e.join(pd.concat([residue_res, group_res], axis = 1 ))

            e = e.fillna(-999)

            e['gap'] = '.'

            try:
                res = np.where(e[0] != -999,
                            np.where(e[0] != '_',
                                       e[0],
                                       np.where((e[0] == '_') & (e[1] != -999),
                                                        e[1],
                                                         e[0]
                                                )
                                     ),
                           np.where(e[1] != -999, e[1], '.')
                           )

            except:
                res = '.' * self.matrix.shape[1]

            dc[threshold] = ''.join(res)


        if output.lower() == 'dict':
            return dc

        elif output.lower() in ['df', 'dataframe', 'data frame']:
            return MsaTable(pd.DataFrame({'consensus_level': [x for x in dc.keys()],
                                          'sequence':[x for x in dc.values()]}),
                            metadata = self.metadata)

        elif output.lower() == 'matrix':
            return self.seq2matrix(df = pd.DataFrame({'consensus_level': [x for x in dc.keys()],
                                                 'sequence':[x for x in dc.values()]},
                                                ), no_self = True)


    def _load_stats(self,sources = [],
                     additional_sources = [],
                     verbose = 0,
                     log_file = '',
                     **kwargs):

        try:
            base = os.path.join(os.path.realpath(os.path.join(os.path.abspath(__file__), '..', '..') ), 'seq')
        except:
            base = os.path.dirname(os.path.realpath(__name__))
        self._verbose = verbose
        self._log_file = log_file
        self.sources =[]
        self.sources.extend(sources)
        self.config = {}
        self._switch_source = loadClasses( base+'.conservation')
        if additional_sources:
            for additional_source in additional_sources:
                self._switch_source.update(loadClasses(additional_source))

    def stats(self, metric = 'shannon_entropy',
                     gap_penalty = 1,
                     pseudocount = 0.0000001,
                     **kwargs
                     ):

        df = self
        self._load_stats()

        if not isinstance(self.matrix, pd.DataFrame):
            self.seq2matrix(**kwargs)

        # Some shortcuts

        if metric == 'se':
            metric = 'shannon_entropy'
        elif metric == 'pe':
            metric = 'property_entropy'
        elif metric == 'jsd':
            metric = 'jensen_shannon_divergence'
        elif metric == 'pre':
            metric = 'property_relative_entropy'
        elif metric == 're':
            metric = 'relative_entropy'

        return self._switch_source[metric](df['sequence'], self.matrix,
                                           gap_penalty = gap_penalty,
                                           pseudocount = pseudocount,
                                           **kwargs).run()

    def conservation_types(self):
        self._load_stats()
        return list(self._switch_source.keys())

    def seq2colors(self):
        df = self
        colors, cmap, bounds, norm = _colors_dict_and_cmap()
        _ = list(df['sequence'].map(lambda x: [int(colors[y])for y in x]))
        return (_, cmap, bounds, norm)

class alignment:
    import warnings
    warnings.filterwarnings('ignore')
    def __init__(self, alignment, ipt_format = 'fasta',
                 **kwargs):

        # Check if alignment is file
        try:
            self.input_file_checker = os.path.isfile(alignment)
        except:
            self.input_file_checker = False
        self.alignment = alignment
        self.ipt_format = ipt_format
        self.kwargs = kwargs

    def read(self, **kwargs):

        if self.ipt_format == 'fasta':
            if self.input_file_checker:
                alignment = open(self.alignment).read().splitlines()
            elif isinstance(self.alignment, str):
                alignment = self.alignment.split('\n')

            else:
                pass

            return MsaTable(self._fasta2df(alignment),
                            metadata = {'sequence': 'sequence',
                                        'Header': 'Header',
                                        'ID': 'ID'})

        elif self.ipt_format in ['table', 'dataframe', 'df']:
            kwargs = self.kwargs
            try:
                metadata = kwargs.pop('metadata', None)

            except:
                pass

            if self.input_file_checker:
                df = pd.read_csv(self.alignment, **self.kwargs)
                if 'sequence' in metadata.keys():
                    df[metadata['sequence']] = df[metadata['sequence']].str.ljust(df[metadata['sequence']].str.len().max(), '-')

                else:
                    df[df.columns[1]] = df[df.columns[1]].str.ljust(df[df.columns[1]].str.len().max(), '-')


                return MsaTable(df,
                                metadata = metadata)
            else:
                return MsaTable(self.alignment, metadata = metadata)

        else:

            idd=[]
            seq=[]

            if (self.input_file_checker):
                fasta_sequences = SeqIO.parse(self.alignment, self.ipt_format)
            else:
                fasta_sequences = SeqIO.parse(iStringIO(self.alignment), self.ipt_format)

            for fasta in fasta_sequences:
                idd.append(fasta.id)
                seq.append(str(fasta.seq))

            data={'ID':idd, 'seq':seq}
            z = pd.DataFrame(data=data)
            # z['len']=z.seq.str.len()

            return MsaTable(z, metadata = {'ID': 'ID',
                                           'sequence': 'seq'})

    def _fasta2df(self, fi, split_seq = False):
        tmp = {}
        splitted = {}
        y = 0
        for x in range(0, len(fi)):
            if fi[x].startswith('>'):
                header = fi[x]
                header = fi[x].replace('>', '')
                tmp[header] = ''
            else:
                tmp[header] += fi[x]
                if split_seq:
                    for e in fi[x]:
                        splitted[str(y)] = e
                        y+=1
        if split_seq:
            df = pd.DataFrame({'Header':list(tmp.keys()),
                               'ID': [x.split(' ')[0] for x in tmp.keys()],
                               'sequence': list(tmp.values())})
            df2 = df['sequence'].apply(lambda x: pd.Series(list(x)))
            df2['Header'] = df['Header']
            df2['ID'] = df['ID']
            return df2
        else:
            df = pd.DataFrame({'Header':list(tmp.keys()),
                               'ID': [x.split(' ')[0] for x in tmp.keys()],
                               'sequence': list(tmp.values())})
            df = df.astype(str)
            # pad
            df['sequence'] = df['sequence'].str.ljust(df['sequence'].str.len().max(), '-')

        return df

def _colors_dict_and_cmap():
    from matplotlib import colors as ccolors
    colors = {'A':1,
                    'I':1,
                    'L':1,
                    'M':1,
                    'F':1,
                    'W':1,
                    'V':1,
                    'K':2,
                    'R':2,
                    'E':3,
                    'D':3,
                    'N':4,
                    'Q':4,
                    'S':4,
                    'T':4,
                    'C':5,
                    'G':6,
                    'P':7,
                    'H':8,
                    'Y':8,
                    'X':0,
                    'a':8, #aromatic
                    'l':1, #alifatic
                    'h':1, #hydrophobic
                    '+':2, #positive
                    'c':8, #charged
                    'p':4 , #polar
                    'o':4, #alcohol
                    'u':1, #tiny
                    's':1, #small
                    'b':1, #big
                    '.':0, #gap
                    '-':3 #negative
                   }

    cmap = ccolors.ListedColormap(['white', (134/256,159/256,234/256), # blue 1
                                                              (220/256,54/256,35/256), # red 2
                                                              (178/256, 80/256,186/256),# magenta 3
                                                              (89/256, 188/256, 60/256), # green 4
                                                              (225/256, 134/256, 131/256), # pink 5
                                                              (227/256, 149/256, 87/256), # orange 6
                                                              (192/256, 192/256, 61/256), # yellow 7
                                                              (76/256, 161/256, 163/256) #cyan 8
                                                              #9 charged
                                                                                        #10 tiny
                                                                                        # 11 small
                               #                                                          # 12 big
                                                                                                                      ])
    bounds=[0,1,2,3,4,5,6,7,8]
    norm = ccolors.BoundaryNorm(bounds, cmap.N)
    return (colors, cmap, bounds, norm)

if __name__ == '__main__':
    print("test read fasta")
    a = alignment('/home/kaihami/projects/ammonium_transp/work/20190208/amt.nr.0d8.amt.domain.sliced.filtered.formated.ordered.modified.trimmed.long_branches_removed.msa').read()


    print('test read dataframe')
    try:
        b = alignment('/home/kaihami/test/fasta/chpT.df', ipt_format = 'df', sep = '\t').read()

    except:
        pass

    b = alignment('/home/kaihami/test/fasta/chpT.df', ipt_format = 'df', header = None,
                 # columns = ['pid', 'seq'],
                  metadata = {},
                  sep = '\t').read()
