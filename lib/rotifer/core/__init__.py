# Set global configuration

import os
import sys
from os.path import expanduser, dirname, realpath

GlobalConfig = {
        'base': dirname(dirname(dirname(dirname(realpath(__file__))))),
        'cache': expanduser(os.path.join(*['~','.cache','rotifer'])),
        'user': expanduser(os.path.join(*['~','.rotifer']))
        }
if not os.path.exists(GlobalConfig['cache']):
    try:
        os.makedirs(GlobalConfig['cache'])
    except:
        print(f'Unable to create cache directory at {GlobalConfig["base"]}', file=sys.stderr)

# Load default modules
from . import cli
from . import functions
from . import loadpath
from . import log
