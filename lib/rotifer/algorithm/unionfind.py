#!/usr/bin/env python3

__version__ = '0.1'

"""
https://en.wikipedia.org/wiki/Disjoint-set_data_structure
"""

## Importing modules

from multiprocessing import Process, Manager
from multiprocessing.managers import BaseManager
import sys
import pandas as pd
import time

## Class
class Disjoint(object):
    """
    Implementation of union-find algorithm
    Peform a union-find over a dataframe
    """
#    def __init__(self):
#        x.p = p
#        x.rank = 0
#    def union(self, x,y):
#        link(find_set(x), find_set(y))
#
#    def link(self, x,y):
#        if x.rank > y.rank:
#            y.p = x
#        else:
#            x.p = y
#            if x.rank == y.rank:
#                y.rank = y.rank +1
#    def find_set(self, x):
#        if x != x.p:
#            x.p = find_set(x.p)
#        return x.p
    def __init__(self,df, colA = '', colB = ''):
        """
        Add dataframe and select columns to perform the union-find
        input data:
            df:   a pandas dataframe
            colA: a source column
            colB: a target column
        Returns two dictionaries
            leader: a dict that could be used in tree search (e.g. Kruskal's algo), contains the element and its root.
            group:  key = tree root, values = components of the tree

        After initializing this class, the user can change the parameter self.sort = [True/False]
        This self.sort will change how the table will be printed (see __str__ method):
            self.sort = True => Ascending order
            self.sort = False => Descending order
        """
        self.tups = list(zip(df[colA], df[colB]))
        self.sort = True
        self.leader = {}
        self.group = {}
        for a,b in self.tups:
           self.add(a,b)

    def add(self, x, y):
        leaderx = self.leader.get(x)
        leadery = self.leader.get(y)
        if leaderx is not None:
            if leadery is not None:
                if leaderx == leadery: return
                groupx = self.group[leaderx]
                groupy = self.group[leadery]
                if len(groupx) < len(groupy):
                    x, leaderx, groupx, y, leadery, groupy = y, leadery, groupy, x, leaderx, groupx
                groupx |= groupy
                del self.group[leadery]
                for k in groupy:
                    self.leader[k] = leaderx
            else:
                self.group[leaderx].add(y)
                self.leader[y] = leaderx
        else:
            if leadery is not None:
                self.group[leadery].add(x)
                self.leader[x] = leadery

            else:
                self.leader[x] = self.leader[y] = x
                self.group[x] = set([x, y])
    def rlist(self):
        '''
        return a list of lists.
        Element of the list correspond to a tree
        example:
            [['A', 'B', 'C'], ['D', 'E']]
            Indicates two trees:
            1: A - B - C
            2: D - E
        '''
        self.ls = []
        for k in self.group.keys():
            _ = list(self.group[k])
            self.ls.append(_)
        return self.ls
    def rdf(self):
        '''
        Returns a pandas dataframe with the following elements:
        Node, id, Components:
        Node is a unique Node name (from the input dataframe)
        id is a internal id, refering to the tree (if two nodes has the same id they belong to the same tree)
        Components is the number of components in the tree
        '''
#        self.rdf = pd.DataFrame(columns = ['Node', 'Root', 'Components'])
        node = []
        root = []
        components = []
        for k in self.group.keys():
            _ = list(self.group[k])
            for ele in _:
                node.append(ele)
                root.append(str(k))
                components.append(len(_))
#                self.rdf.loc[self.rdf.shape[0]] = [ele, str(k), len(_)]
        dc = {'Node': node,
              'Root': root,
              'Components': components}
        self.rdf = pd.DataFrame.from_dict(dc) #columns = ['Node', 'Root', 'Components'])
        dsize = pd.DataFrame({'count': self.rdf.groupby('Root').size()}).sort_values(by='count', ascending=False).reset_index()
        dsize['id'] = dsize.index
        m = dict(zip(dsize['Root'], dsize['id']))
        self.rdf['id'] = self.rdf['Root'].map(m)
        self.rdf = self.rdf.drop(columns = ['Root'] )
        self.rdf.sort_values(['Components'])
        cols = ['Node', 'id', 'Components']
        self.rdf = self.rdf[cols]
        self.rdf.sort_values(['Components'])
        self.rdf = self.rdf.sort_values(['Components'])
        return self.rdf

    def __str__(self):
        """
        usage: print(object)
        Print the union find table separated by tabs (\\t)
        Node, id, Components
        From the lowest to highest number of components
        Reverse, set object.sort = False before calling print
        example:
        > uf = Disjoint(df,'colA','colB')
        > print(uf)
        >> WP_017899391.1  2       285
        >> WP_043081594.1  2       285
        >> WP_061066926.1  2       285
        >> WP_003017908.1  2       285
        >> WP_001556401.1  1       300
        >> WP_001366798.1  1       300
        >> WP_024130714.1  1       300

        > uf.sort = False
        > print(uf)
        >> WP_001556401.1  1       300
        >> WP_001366798.1  1       300
        >> WP_024130714.1  1       300
        >> WP_017899391.1  2       285
        >> WP_043081594.1  2       285
        >> WP_061066926.1  2       285
        >> WP_003017908.1  2       285
        """
        node = []
        root = []
        components = []
        for k in self.group.keys():
            _ = list(self.group[k])
            for ele in _:
                node.append(ele)
                root.append(str(k))
                components.append(len(_))
#                self.rdf.loc[self.rdf.shape[0]] = [ele, str(k), len(_)]
        dc = {'Node': node,
              'Root': root,
              'Components': components}
        self.rdf = pd.DataFrame.from_dict(dc) #, columns = ['Node', 'Root', 'Components'])
        dsize = pd.DataFrame({'count': self.rdf.groupby('Root').size()}).sort_values(by='count', ascending=False).reset_index()
        dsize['id'] = dsize.index
        m = dict(zip(dsize['Root'], dsize['id']))
        self.rdf['id'] = self.rdf['Root'].map(m)
        self.rdf = self.rdf.drop(columns = ['Root'] )
        self.rdf.sort_values(['Components'])
        cols = ['Node', 'id', 'Components']
        self.rdf = self.rdf[cols]
        self.rdf = self.rdf.sort_values(['Components'], ascending = self.sort)
        return (self.rdf.to_csv(sep = "\t", index = None))

