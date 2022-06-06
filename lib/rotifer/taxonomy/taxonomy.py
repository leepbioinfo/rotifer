import pandas as pd
import warnings
import numpy as np
warnings.filterwarnings('ignore')
class df2tax:
    '''
    Classification made easy.

    Usage:
    > import sys
    > sys.path.insert(0,'/home/kaihami/mymodules/')
    > from rotifer.taxonomy import taxonomy
    > taxonomy.df2tax()
    '''
    def __init__(self, df, column = 'classification'):
        '''
        Load a rneighbors dataframe
        Parameters:
        df: A pandas dataframe with a classification column
        column: the column name with NCBI lineage
        '''
        self.df = df
        self.annotation= pd.read_csv('/home/kaihami/projects/gen_scripts/ncbi_tax/taxdump.ncbi.lineage.tsv', sep = '\t')

        self.column = column
        self.taxonomy = self._prepare(self.column)

    def _prepare(self, column = 'classification'):
        self.subdf = pd.DataFrame({'name': np.concatenate(self.df[column].str.split(';').str[:-1].values),
                              'idx': self.df.index.repeat(self.df[column].str.split(';').str[:-1].apply(len))
                             })

        self.subdf['name'] = self.subdf['name'].map(lambda x: x.strip())

        finaldf = self.subdf.merge(self.annotation, right_on = 'name_txt', left_on = 'name', how = 'left')
        _2 = finaldf[finaldf['rank'].isin(['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])]
        checker = {'superkingdom': 0,
                   'phylum': 1,
                    'class': 2,
                     'order': 3,
                     'family': 4,
                  'genus': 5,
                  'species': 6
                  }
        g = _2.groupby('idx').size()

        _2 = _2.drop(_2.index[g[g==1].index])
        _2['format'] = np.where((_2.shift(-1)['name'] == _2['name']) | (_2['name'] == _2.shift()['name']) , 'eq','neq')
        _2['format2'] = np.where((_2['format'] == 'eq'),
                                np.where(_2.shift()['format'] != 'eq', _2.shift()['rank'], _2.shift(2)['rank']), 'no')
        sub = _2[_2['format'] == 'eq']
        sub['current'] = sub['rank'].map(checker)
        sub['test'] = sub['format2'].map(checker)
        sub['dif'] = sub['current'] - sub['test']
        remove_idx = sub[sub['dif']!= 1].index

        finaldf = finaldf.drop(remove_idx)

        finaldf.drop_duplicates(inplace = True)
        finaldf= finaldf[['name', 'idx','rank']]

        return finaldf

    def class_level(self):
        '''
        Show possibles taxomony levels
        '''
        print('\n'.join([str(x) for x in list(self.taxonomy['rank'].unique())]))

    def class_list(self):
        '''
        Show possibles tax level list format
        '''
        return ([str(x) for x in list(self.taxonomy['rank'].unique())])

    def tax(self, tax = ''):
        '''
        Pass the taxonomy level (e.g. phylum)
        '''

        result = self.taxonomy[self.taxonomy['rank'] == tax]
        return result[['idx','name', 'rank']]

    def finder_name(self, name = ''):
        '''
        Find using classification name (e.g. Alphaproteobacteria)
        '''


        result = self.taxonomy[self.taxonomy['name'] == name]
        #result = result.drop_duplicates()
        return result[['idx','name', 'rank']]

    def find_regex(self, regex = '', case_sensitive = True):
        '''
        Use regex to find ocurrence (case_sensitive = True/False, if False bacteria/Bacteria results in the same output)
        '''
        if case_sensitive:
            idx = self.taxonomy['name'].str.contains(regex, regex=True)
        else:
            idx = self.taxonomy['name'].str.contains(regex, regex=True, flags=re.IGNORECASE)

        result = self.taxonomy[idx]

        return result[['idx','name', 'rank']]

    def rlineage(self, levels = ['superkingdom', 'phylum','class']):
        '''
        Input a list of taxonomic levels and return in our favourite format.
        INPUT:
        levels: TAXONOMIC LEVELS [LIST]

        OUTPUT:
        A dataframe containg indexes and lineages separated by '>'

        Example:
        > taxonomy = df2tax(ipt)
        > taxonomy.rlineage()
        ... idx    lineage
        ...   0    Bacteria>Acidobacteria>Acidobacteriia
        '''

        tmp = self.taxonomy[self.taxonomy['rank'].isin(levels)]
        tmp['name'] = self.taxonomy.name + '>'
        g = tmp.groupby('idx').name.sum()
        g = g.to_frame()
        g.columns = ['lineage']
        g['lineage'] = g['lineage'].str.rstrip('>')

        return g
