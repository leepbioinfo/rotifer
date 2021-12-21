# This module's documentation
__doc__ = '''Rotifer's core configuration and setup'''

# Dependencies
import os
import sys

# Global configuration
def __setup_base_config():
    from os.path import expanduser, dirname, realpath

    # Base directories
    baseDir = dirname(dirname(dirname(dirname(realpath(__file__)))))
    userDir = expanduser(os.path.join(*['~','.rotifer']))

    # Main dictionary: this is the most important part of this module!!!!
    GlobalConfig = {
            'base': baseDir,
            'baseConfig': os.path.join(baseDir,"etc","rotifer"),
            'cache': expanduser(os.path.join(*['~','.cache','rotifer'])),
            'user': userDir,
            'userConfig': os.path.join(userDir,"etc"),
            }

    # Create directiories, if needed
    auto = [ 'cache', 'user', 'userConfig' ]
    for target in sorted([ GlobalConfig[x] for x in auto ]):
        if not os.path.exists(target):
            try:
                os.makedirs(target)
            except:
                print(f'Unable to create cache directory at {target}', file=sys.stderr)

    # Return main configuration
    return GlobalConfig

# Setup configuration
GlobalConfig = __setup_base_config()

# Load default modules
from . import cli
from . import functions
from . import loadpath
#from . import log

