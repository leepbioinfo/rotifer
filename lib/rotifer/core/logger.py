import os
import sys
import logging
from tqdm import tqdm

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
def get_logger(name):
    RotiferDefaultFormatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
    RotiferDefaultHandler = TqdmLoggingHandler()
    RotiferDefaultHandler.setFormatter(RotiferDefaultFormatter)
    logger = logging.getLogger(name)
    logger.addHandler(RotiferDefaultHandler)
    return logger

