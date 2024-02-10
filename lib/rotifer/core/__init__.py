# This module's documentation
__doc__ = '''Rotifer's core configuration and setup'''

# Dependencies
import os
import sys

# Core configuration
def _get_core_config():
    from os.path import expanduser, dirname, realpath

    # Base directories
    baseDir = dirname(dirname(dirname(dirname(realpath(__file__)))))
    userDir = expanduser(os.path.join('~','.rotifer'))

    # Main dictionary: this is the most important part of this module!!!!
    config = {
        'base': baseDir,
        'baseConfig': os.path.join(baseDir,"etc","rotifer"),
        'baseDataDirectory': os.path.join(baseDir,"share","rotifer"),
        'user': userDir,
        'userConfig': os.path.join(userDir,"etc"),
        'userDataDirectory': os.path.join(userDir,"share"),
    }

    # Libraries
    libdir = os.path.join(userDir,"lib","python")
    if libdir not in sys.path:
        sys.path.insert(0,libdir)

    return config

config = _get_core_config()
