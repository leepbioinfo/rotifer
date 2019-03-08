#!/usr/bin/env python3

__version__ = "0.01"

### Import rotifer package
import sys
import os
sys.path.insert(0, os.path.join(os.getcwd(), '../..'))

### Import core cli
import rotifer.core.cli as corecli

### Other packages
import argparse
import sys

### create arguments for sequences


class action:
    '''
    Rparser custom actions!
    '''
    def openload(self, parser, largs):
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error('No input!')
        elif not '-' in largs:
            largs.insert(0, sys.stdin)
        sys.stdin


    class autoopen_seqio(argparse.Action):
        '''
        An action to decide file type, open, and load the file into memory, close after
        So far it decide if string, file, stdin
        Optional parameters:
            duplicates = Boolean [True/False]
        '''
        def __init__(self, option_strings, file_type = 'list', duplicates = True,
                      *args, **kwargs):
            self._add1 = file_type
            self._add2 = duplicates
            super(action.autoload, self).__init__(option_strings=option_strings,
                                           *args, **kwargs)
        def __call__(self, parser, namespace, values, option_string = None):
            #TODO insert new function for the var file type
            #print(self._add1)
            if self._add2 == True:
                f = action.openread(self,parser, values)
            if self._add2 == False:
                f = action.openread2(self,parser, values)
            setattr(namespace, self.dest, f)

### Sequence actions




