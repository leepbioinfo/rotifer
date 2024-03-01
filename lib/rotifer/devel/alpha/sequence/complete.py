#!/usr/bin/env python3

def complete(df, include=True, values=False):
    import pandas as pd
    import rotifer.interval.utils as riu
    complete = pd.concat([df.copy(),riu.complement(df.copy(), include, values)]).sort_values(['ID','start','end'])
    return complete


