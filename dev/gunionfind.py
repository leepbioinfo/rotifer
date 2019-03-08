import os
import sys
sys.path.insert(0, os.path.join('/home/kaihami/mymodules'))
import rotifer.algorithm.unionfind as union
#import rotifer.algorithm.pure.unionfind as union

import pandas as pd
import time
if __name__ == '__main__':
    df = pd.read_csv(sys.argv[1], sep = '\t', header = None)
#    print(df.head())
#    df.columns = ['a', 'b']
#    print(df['a'])
#    print(df)
    s = time.time()
    uf = union.Disjoint(df, 0,1)
    uf.sort = False
    print(uf)
    e = time.time()
    print('')
    print('#'*20)
    print(e-s)
    #uf.rprint()
