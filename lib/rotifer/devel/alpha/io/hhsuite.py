# HH-suite parser

__author__ = "Robson Francisco de Souza, Rodolfo Alvarenga Ribeiro"
__version__ = 0.001

# Libraries
import os
import sys
import pandas as pd
from glob import glob
from rotifer.core.functions import import_path

# Add path to rotifer libraries to Python's library search list
hhr = import_path(os.path.realpath(__file__).split("rotifer")[0] + "rotifer/bin/hhsuite2table")

# Parse output tables
def read_hhr(indir, each=False):
    '''
    Parse all HH-suite (*.hhr) output files located in a
    directory or, alternatively, a list of .hhr files

    Usage:
      # Parsing all files at once
      from rotifer.devel.alpha.io import read_hhr
      df = read_hhr(hhrdir)

      # Iterating over each file
      from rotifer.devel.alpha.io import read_hhr
      for hhr in read_hhr(hhrdir):
        hhr.to_csv(hhr.replace("hhr","tsv"), sep="\t", index=False)

    Arguments:
      - indir : a directory with .hhr files or a list of .hhr files
      - each  : instead of parsing all files at once, setting this
                options to true will create a generator that will
                parse and return the contents of each .hhr file at
                each iteration

    Returns:
      A Pandas dataframe or a generator of Pandas dataframes
    '''
    if not isinstance(indir,list):
        if os.path.exists(indir) and os.path.isdir(indir):
            indir = os.path.realpath(indir)
            indir = glob(f'{indir}/*.hhr')
        else:
            print(f'{__name__}.read_hhr: unknown directory {indir}', file=sys.stderr)
            return None
    if each:
        return __parse_each(indir)
    else:
        return __parse_all(indir)
    return df

def __parse_each(indir):
    for infile in indir:
        yield hhr.hhsuite2pandas(infile)

def __parse_all(indir):
    return hhr.hhsuite2pandas(indir)
