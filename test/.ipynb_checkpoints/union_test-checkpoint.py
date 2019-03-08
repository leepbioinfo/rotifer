#/usr/bin/env python3

import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.algorithm.unionfind2 as union
#import rotifer.algorithm.pure.unionfind as union

import pandas as pd
import time
if __name__ == '__main__':
    df = pd.read_csv('./data/union_find_df.tsv', sep = '\t', header = None)
    s = time.time()
    # create an empty graph
    uf = union.Disjoint()
    load_list = [['A','B'],
                 ['B', 'C'],
                 ['A', 'C'],
                 ['D', 'E']]
    uf.from_list(load_list)
    a = uf.add('G','H')
    ### Return from load_list => A,B,C ; D,E
#   Node    id      Components
#   B       0       3
#   C       0       3
#   A       0       3
#   E       1       2
#   D       1       2
    uf.sort = True
    print('## List test')
    print(uf)
    uf2 = union.Disjoint()
    uf2.from_df(df, 0,1)
    uf2.sort = False
    print('## DataFrame no weight')
    print(uf2.rdf().tail())

    uf3 = union.Disjoint()
    uf3.from_df_weight(df, 0,1,2, higher = 0.5, lower = 0.5)
    uf3.sort = False
    print('## DataFrame with weight')
    print(uf3.rdf().tail())


#    e = time.time()
#    print('')
#    print('#'*20)
#    print(e-s)
    #uf.rprint()

