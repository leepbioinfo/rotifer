#!/usr/bin/env python3
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import pandas as pd
import numpy as np
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from multiprocessing import Pool, Process
import itertools
import shutil
class cluster_neighbor:
    def __init__(self, df, blockid = '', cluster = ''):
        '''
        You can access the data count using .matrix
        '''
        self.gb = df.groupby(blockid)[cluster].apply(lambda x: x.tolist()).to_dict()
        self.threads = 3
        self.matrix = pd.DataFrame(columns = list(self.gb.keys()))
#                                   columns = list(self.gb.keys()))
        size = int(len(self.gb.keys())/self.threads)
        sub_chunck  = [list(self.gb.keys())[i:i+size] for i in range(0, len(self.gb.keys()), size)]
        self.chunk = {index: x for index, x in enumerate(sub_chunck, start=1)}
        _ = [sorted(p) for p in itertools.product(list(self.chunk.keys()), repeat=2)]
        self.combinations = []
        for ele in _:
            if ele not in self.combinations:
                self.combinations.append(ele)
            else: pass

        print(self.combinations)
    def _sub_matrix(self, args):
        '''
        Build a square matrix
        '''
        vals = args
        print(vals)
        val1, val2 = self.chunk[vals[0]], self.chunk[vals[1]] # two lists
        name = 'tmp_matrix_'+str(vals[0])+'_'+str(vals[1])
        tmp_matrix = pd.DataFrame(index = val1, columns = val2)
        for i in val1:
            for j in val2:
                # compare
                array_i = self.gb[i]
                array_j = self.gb[j]
                res = len([x for x in array_i if x in array_j])
                tmp_matrix.loc[i,j] = res
        tmp_matrix.to_csv('./tmp/'+name, sep = '\t')
    def _combine_matrix(self):
        for _matrix in os.listdir('./tmp'):
            _tmpload = pd.read_csv('./tmp/'+_matrix, sep = '\t', index_col = 0)
            _tmpmap =  _tmpload.groupby(blockid)[cluster].apply(lambda x: x.tolist()).to_dict()
            print(_tmpload.head())
            self.matrix = self.matrix.join(_tmpload, how='outer', on = _tmpload.columns.tolist())
            print(self.matrix)
        return self.matrix
    def clustering(self):
        '''
        Get linkage
        '''
        correlation_array = np.asarray(self.matrix)
        self.linkage = hc.linkage(correlation_array, method = 'average')
#        self.linkage = hc.linkage(sp.distance.squareform(correlation_array), method = 'average')
        return self.linkage
    def plt(self):
        '''
        Plot
        '''
        sns_plot = sns.clustermap(self.matrix, row_linkage=self.linkage, col_linkage=self.linkage)
        sns_plot.savefig('aa.png')

    def build_matrix(self):
        try:
            print(os.getcwd())
            os.mkdir('tmp')
        except:
            print('Already tmp')
        self.number_process = 1
        jobs = []
        for x in range(len(self.combinations)):
            self.p = Process(target = self._sub_matrix, args = (self.combinations[x],))
            jobs.append(self.p)
            self.p.start()

        for p in jobs:
            p.join()
        self._combine_matrix()
        try:
            shutil.rmtree('./tmp/')
        except:
            print('Not removed')

