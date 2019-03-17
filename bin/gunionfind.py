#!/usr/bin/env python3

# Rotifer's imports
import os
import sys
_d = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
for _d in [ os.path.join(_d, "lib"), os.path.join(_d, "lib", "python" + str(sys.version_info.major) + "." + str(sys.version_info.minor), "site-packages") ]:
    if os.path.exists(_d):
        sys.path.insert(0,_d)
import rotifer.core.cli as corecli
import rotifer.algorithm.unionfind as union
#import rotifer.algorithm.pure.unionfind as union

# Other imports
import pandas as pd
import time

# Main program
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
