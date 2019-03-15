import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))


import rotifer.core.cli as cli
print(cli.my())
'''
from rotifer.core import cli

print((rotifer.core.cli.my()))
'''
