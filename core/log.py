#!/usr/bin/env python3

from datetime import datetime as dt
import mdv
from os.path import expanduser
import yaml
import os
import sys

__version__ = 0.1
__authors__ = 'Gilberto Kaihami; Robson Souza'
import socket

def log(message = {}, level = 0, log_file = '', name = ''):
    """
    Output a message to stderr.

    Three options available:
    1-) Input a string
    2-) Input a list with one element
    3-) Input a list with n elements

    If a list with more than one element is passed,
    the first element will be printed with timelog the other lines will be printed in a new line.
    -----------
    PARAMETERS:
    message:  A dictionary containing the messages.
    level:    Verbose level (0-3) (default: 0)
    log_file: A file to save the log (default: '')
    name:     Script name (default: __name__)
    -----------
    OUTPUT:
    A message in the format for option (1) and (2)
    ## [DAY/MONTH Hour:Minute:Second] <message>
    or

    A message in the format for option (3)
    ## [DAY/MONTH Hour:Minute:Second] <message>
    ### <message>
    ### ...
    -----------
    Example:
    from rotifer.core.log import log

    if verbose:
        log({1: 'Info message',
             2: 'Warning message',
             3: 'Debug message'},
             level = verbose,
             log_file = log_file,
             )

    if verbose:
        log({1: 'Info message',
             2: 'Warning message',
             3: 'Debug message'},
             level = verbose,
             log_file = log_file,
             name = __name__ # only needed if used inside a package
             )

    """
    import __main__ as mainprog
    import sys
    now = dt.now().strftime('%D %H:%M:%S')
    max_level = 3
    hostname = socket.gethostname()

    _debug_switch = {1: 'INFO',
                     2: 'WARNING',
                     3: 'DEBUG'}

    if level > max_level:
        level = max_level
    else:
        pass
    if name == '__main__':
        name = (sys.argv[0]).replace('.py', '')
    else:
        name = (os.path.basename(name))


    col_len = _collen(now, hostname, name)
    try:
        for x in range(1,level+1):

            try:
                if isinstance(message[level], list):
                    msg_print = ' '.join([x.strip() for x in message[x]])

                else:
                    msg_print = message[x].strip()

                s = [now, hostname, name, _debug_switch[x],':', msg_print]
                if not log_file:
                    to_print = [s[int(i)].ljust(col_len[int(i)]) for i in range(0, len(s))]
                    sys.stderr.write('\t'.join(to_print) + '\n')
                    sys.stderr.flush()
                else:
                    with open(log_file, 'a') as f:
                        to_print = [s[int(i)].ljust(col_len[int(i)]) for i in range(0, len(s))]
                        f.write('\t'.join(to_print) + '\n')
            except:
                pass

    except:
        pass

def _collen(now, hostname, name, levels = ['INFO', 'WARNING', 'DEBUG']):
    collen_len = []
    collen_len.append(len(now))
    collen_len.append(len(hostname))
    collen_len.append(len(name))
    collen_len.append(max([len(x) for x in levels]))
    collen_len.append(0)
    collen_len.append(0)
    return collen_len

if __name__ == '__main__':
    import t
    log({1:'This is info'}, level = 1)
    print()
    log({1: 'This info',
          2: 'This is Warning',
          3: 'DEBUG'},
         level = 1)

    print()
    log({1: 'This info',
          2: 'This is Warning',
          3: 'DEBUG'},
         level = 2)
    #
    # # Exceed
    print()
    log({
          3: 'DEBUG 3'},
         level = 3)
    print()
    log({1: 'This info',
          2: 'This is Warning',
          3: 'DEBUG exce'},
         level = 4)
    print()
    #
    # # Write to log file
    log({1: 'This info',
          2: 'This is Warning',
          3: 'DEBUG'},
         level = 3,
         log_file = 'a')