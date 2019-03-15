#!/usr/bin/env python3

# For now, we use this as part of our script's template
import os
import sys
_add_path = [ os.path.join(os.path.dirname(os.path.dirname(__file__)), "lib"),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "lib", "python" + str(sys.version_info.major) + "." + str(sys.version_info.minor), "site-packages")
        ]
for _d in _add_path:
    if os.path.exists(_d):
        sys.path.insert(0,_d)

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
