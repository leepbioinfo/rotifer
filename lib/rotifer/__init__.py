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
import logging as logging
from rotifer.core import logger as rcl

def _setup_config(logger):
    from rotifer.core import config as CoreConfig 
    from rotifer.core.functions import loadConfig

    # Reset configuration using environment variables
    config = loadConfig(":" + __name__, defaults = {
        'data': os.environ['ROTIFER_DATA'] if 'ROTIFER_DATA' in os.environ else "/databases",
        **CoreConfig
    })

    # Create directiories, if needed
    auto = [ 'cache', 'user', 'userConfig' ]
    for target in sorted([ config[x] for x in auto ]):
        if not os.path.exists(target):
            logger.warning(f'Creating directory {target} ...')
            try:
                os.makedirs(target)
            except:
                logger.warning(f'Unable to create cache directory at {target}')

    # Return main configuration
    if 'cache' not in config:
	config['cache'] = '/tmp'
    if not os.path.exists(config['cache']):
        config['cache'] = os.path.join(os.path.expanduser("~"),".cache","rotifer")
    return config

# Setup configuration
logger = rcl.get_logger(__name__)
config = _setup_config(logger)
GlobalConfig = config
