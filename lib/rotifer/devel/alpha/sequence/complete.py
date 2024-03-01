#!/usr/bin/env python3

import rotifer
from rotifer import GlobalConfig
logger = rotifer.logging.getLogger(__name__)

def complete(df, include=True, exclude=False, values=False):
    import pandas as pd
    import rotifer.interval.utils as riu
    complete = pd.concat([df, riu.complement(df, include, exclude, values)]).sort_values(['ID','start','end'])
    return complete


