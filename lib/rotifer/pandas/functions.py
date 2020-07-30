#!/usr/bin/env python3

import os
import sys
import logging
import pandas as pd

def print_everything(max_rows=10000000000, max_columns=1000000000, max_colwidth=1000000000, width=100000, verbose=True):
    """
    A small function to allow pandas to print not only a few columns and rows.
    Take care when using it: blarge dataframes could take a long time to print.
    """
    pd.options.display.max_rows = max_rows
    pd.options.display.max_columns = max_columns
    pd.options.display.max_colwidth = max_colwidth
    pd.options.display.width = width
    if verbose:
        print('Pandas printing the full dataframe')
