__doc__ = '''
Rapid Open-source Tools and Infrastructure For data Exploration and Research
============================================================================

Rotifer is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.
'''

# Import modules
from tqdm import tqdm
import logging as logging
from rotifer.core import GlobalConfig
from rotifer.core.functions import loadConfig

# Sett Rotifer's root logger

# Classes
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
RotiferDefaultFormatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
RotiferDefaultHandler = TqdmLoggingHandler()
RotiferDefaultHandler.setFormatter(RotiferDefaultFormatter)
logger = logging.getLogger(__name__)
logger.addHandler(RotiferDefaultHandler)
