#!/usr/bin/env python3

import pandas as pd

class test(pd.DataFrame):
    t = 1
    _metadata = ['t']

    def __init__(self,t =2, *args, **kwargs):
        super(test, self).__init__(*args, **kwargs)
        self.t = 2

    @property
    def _constructor(self):
        return test
