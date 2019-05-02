#!/usr/bin/env python3

import os
import warnings
warnings.filterwarnings("ignore")
import sys
from subprocess import PIPE, Popen
from prettytable import PrettyTable

def blockPrint():
    sys.stderr = open(os.devnull, 'w')
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

if os.path.islink(__file__):
    folder = os.path.dirname(os.path.realpath(__file__))

else:
    folder = os.path.dirname(os.path.abspath(__file__))

t = PrettyTable(['Name', 'Description'])

res = []
for p in os.listdir(folder):
    path = os.path.join(folder, p)
    if os.path.isfile(path) and 'rprograms' not in p:
        try:
            os.chdir(folder)
            try:
                blockPrint()
                sub = __import__(p.replace('.py', ''), fromlist = ['__rdoc__'])
            except:
                pass
            description = sub.__rdoc__.split('\n')
            s = ''

            for x in range(0, len(description)):
                if 'DESCRIPTION' in description[x]:
                    for y in range(x+1, len(description)):
                        if description[y].isupper():
                            break
                        else:
                            s += description[y]
                    break
            try:
                res.append([p.replace('.py', ''), s])
            except:
                pass

        except:
            pass
    else:
        pass

enablePrint()
for e in res:
    t.add_row([e[0], e[1]])
print(t)
