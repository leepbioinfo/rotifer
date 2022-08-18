__doc__ = '''
Rapid Open-source Tools and Infrastructure For data Exploration and Research
============================================================================

Rotifer is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.
'''

# Import modules
import os
import sys
from tqdm import tqdm
import logging as logging

# Rotifer's root logger handler class
# https://stackoverflow.com/questions/38543506/change-logging-print-function-to-tqdm-write-so-logging-doesnt-interfere-wit
class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

# Create root logger and default handlers
def _logger():
    RotiferDefaultFormatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
    RotiferDefaultHandler = TqdmLoggingHandler()
    RotiferDefaultHandler.setFormatter(RotiferDefaultFormatter)
    logger = logging.getLogger(__name__)
    logger.addHandler(RotiferDefaultHandler)
    return logger

# Global configuration
def _setup_base_config(logger):
    from os.path import expanduser, dirname, realpath

    # Base directories
    baseDir = dirname(dirname(dirname(realpath(__file__))))
    userDir = expanduser(os.path.join('~','.rotifer'))

    # Main dictionary: this is the most important part of this module!!!!
    GlobalConfig = {
        'base': baseDir,
        'baseConfig': os.path.join(baseDir,"etc","rotifer"),
        'baseDataDirectory': os.path.join(baseDir,"share","rotifer"),
        'cache': expanduser(os.path.join('~','.cache','rotifer')),
        'data': os.path.join(baseDir,"share","rotifer"),
        'user': userDir,
        'userConfig': os.path.join(userDir,"etc"),
        'userDataDirectory': os.path.join(userDir,"share"),
    }

    # Create directiories, if needed
    auto = [ 'cache', 'user', 'userConfig' ]
    for target in sorted([ GlobalConfig[x] for x in auto ]):
        if not os.path.exists(target):
            try:
                os.makedirs(target)
            except:
                logger.error(f'Unable to create cache directory at {target}', file=sys.stderr)

    # Reset configuration using environment variables
    envmap = {'ROTIFER_DATA':'data'}
    for name, value in envmap.items():
        if name in os.environ:
            GlobalConfig[value] = os.environ[name]

    # Return main configuration
    return GlobalConfig

# Setup configuration
logger = _logger()
GlobalConfig = _setup_base_config(logger)
