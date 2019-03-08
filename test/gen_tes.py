import random

_ = open('a').read().splitlines()

for line in _:
    a = line.split('\t')[0:2]
    r = random.uniform(0,1)
    print('\t'.join(a+[str(r)]))
