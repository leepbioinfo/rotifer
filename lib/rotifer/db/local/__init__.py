__doc__ = """
Rotifer connections to local databases
======================================
"""

import re
import os
import types
import typing
import sqlite3
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import StringIO

import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
from rotifer.db.core import BaseCursor
from rotifer.genome.data import NeighborhoodDF
from rotifer.genome.utils import seqrecords_to_dataframe
import rotifer.devel.beta.sequence as rdbs
logger = rotifer.logging.getLogger(__name__)

# Defaults
config = loadConfig(__name__.replace('rotifer.',':'), defaults = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
})
