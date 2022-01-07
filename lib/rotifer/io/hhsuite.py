# HH-suite parser

__author__ = "Robson Francisco de Souza, Rodolfo Alvarenga Ribeiro"
__version__ = 0.002

# Libraries
import os
import sys
import pandas as pd
from glob import glob
from rotifer.core.functions import import_path

# Add path to rotifer libraries to Python's library search list
hhr = import_path(os.path.realpath(__file__).split("rotifer")[0] + "rotifer/bin/hhsuite2table")

# Parse output tables
def read_hhr(indir):
    '''
    Parse **all** HH-suite (*.hhr) output files located in a
    directory or, alternatively, a list of .hhr files

    Parameters
    ----------
    indir : a directory with .hhr files or a list of .hhr files
    
    Returns
    -------
      A Pandas dataframe
    
    Examples
    --------
    Parsing all files at once

    >>> from rotifer.devel.alpha.io import read_hhr
    >>> df = read_hhr(hhrdir)
    '''
    if not isinstance(indir,list):
        if os.path.exists(indir) and os.path.isdir(indir):
            indir = os.path.realpath(indir)
            indir = glob(f'{indir}/*.hhr')
    return hhr.hhsuite2pandas(indir)

def parse_hhr(indir):
    '''
    Parse HH-suite (*.hhr) output files.

    Parameters
    ----------
    indir : a directory with .hhr files or a list of .hhr files

    Returns
    -------
      A generator of Pandas dataframes

    Examples
    --------
    Iterating over each file

    >>> from rotifer.devel.alpha.io import parse_hhr
    >>> for hhr in parse_hhr(hhrdir):
    >>>    hhr.to_csv(hhr.replace("hhr","tsv"), sep="\t", index=False)
    '''
    if not isinstance(indir,list):
        if os.path.exists(indir) and os.path.isdir(indir):
            indir = os.path.realpath(indir)
            indir = glob(f'{indir}/*.hhr')
    for infile in indir:
        yield hhr.hhsuite2pandas(infile)
