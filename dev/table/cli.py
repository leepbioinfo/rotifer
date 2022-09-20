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

### Create arguments for table
class table:
    def __init__(self):
        self.header = ['y', 'header', '']
        self.include = ['i', 'include', '']
        self.exclude = ['e', 'exclude', '']
        self.filter = ['f', 'filter', '']
        self.order = ['o', 'order', '']

        self.helpdc = {'header': '',
                       'include': '',
                       'exclude': '',
                       'filter': '',
                       'order': ''}

    def input(self, number = 1, exclude_filter = [], helper = {}):
        '''
        Set all parameters for a input table
        Customize:
            - exclude_filter = A list with parameters to exclude from the table
            - helper: Add custom help message for the key
        '''
        self.parser = argparse.ArgumentParser(add_help = False)

        for k,v in self.helpdc.items():
            if k in helper.keys():
                self.helpdc[k] = helper[k]


        parameters = {'header': self.header,
                      'include': self.include,
                      'exclude': self.exclude,
                      'filter': self.filter}
        for k, v in parameters.items():
            parameters[k] = ['-'+v[0], '--'+v[1], self.helpdc[k] ]

        if number == 1:
            for k, v in parameters.items():
                if k not in exclude_filter:
                    self.parser.add_argument(v[0],v[1], help = v[2])

        else:
            for v in parameters.values():
                self.parser.add_argument(v[9], v[1], help = v[2])
            for x in range(1, number +1):
                for k, v in parameters.items():
                    shortName = v[0]+str(x)
                    longName = v[1]+str(x)
                    self.parser.add_argument(shortName, longName, help = v[2])

        return self.parser

    def output(self, number = 1, exclude_filter = []):
        self.parser = argparse.ArgumentParser(add_help = False)
        parameters = {'header': self.header,
                      'include': self.include,
                      'exclude': self.exclude,
                      'filter': self.filter}
        for k,v in parameters.items():
            parameters[k] = ['-o'+v[0], '--output_'+v[1], v[2]]
        if number == 1:
            for k, v in parameters.items():
                if k not in exclude_filter:
                    self.parser.add_argument(v[0],v[1], help = v[2])
        return self.parser

def parseArgs(*args):
    '''
    Put all arguments together
    '''
    # print(args)
    if len(args) == 1:
        return argparse.ArgumentParser(parents=args).parse_args()
    else:
        s = argparse.ArgumentParser(parents=args, add_help = True)
        return s.parse_args()

class action:
    '''
    Rparser custom actions!
    '''
    def openread(self,parser,largs):
        '''
        Open, read, and close a file. Remove duplicates
        '''
        # Process sys.stdin
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error("No input!")
        elif not "-" in largs:
            largs.insert(0,sys.stdin)

        # Process sys.argv
        myset=set()
        for iohandle in largs:
            try:
                iohandle = open(iohandle, mode="r")
                myset=myset.union(iohandle.read().splitlines())
                iohandle.close
            except:
                try:
                    myset=myset.union(iohandle.read().splitlines())
                except:
                    if iohandle == "-":
                        if sys.stdin.isatty():
                            parser.error("No input!")
                        else:
                            myset=myset.union(sys.stdin.read().splitlines())
                    else:
                        myset.add(iohandle)
        return(myset)

    def openread2(self, parser, largs):
        '''
        Open, read, and close a file. Keep duplicates.
        '''
        # Process sys.stdin
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error("No input!")
        elif not "-" in largs:
            largs.insert(0,sys.stdin)

        # Process sys.argv
        myset= []
        for iohandle in largs:
            try:
                iohandle = open(iohandle, mode="r")
                myset.extend(iohandle.read().splitlines())
                iohandle.close
            except:
                try:
                    myset.extend(iohandle.read().splitlines())
                except:
                    if iohandle == "-":
                        if sys.stdin.isatty():
                            parser.error("No input!")
                        else:
                            myset.extend(sys.stdin.read().splitlines())
                    else:
                        myset.extend(iohandle)
        return(myset)

    def openload(self, parser, largs):
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error("No input!")
        elif not "-" in largs:
            largs.insert(0,sys.stdin)

        # Process sys.argv
        mylist = list()
        for iohandle in largs:
            try:
                iohandle = open(iohandle, mode="r")
                mylist.append(iohandle)
            except:
                try:
                    mylist.append(iohandle)
                except:
                    if iohandle == "-":
                        if sys.stdin.isatty():
                            parser.error("No input!")
                        else:
                            mylist.append(iohandle)
                    else:
                        mylist.append(iohandle)
        return mylist

    class autoload(argparse.Action):
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

    class autoopen(argparse.Action):
        '''
        An action to decide file type, and open.
        '''
        def __call__(self, parser, namespace, values, option_string = None):
            f = action.openload(self, parser, values)
            setattr(namespace, self.dest, f)

class add_elements:
    def __init__(self):
        self.dc = {}
        self.name = 1
    def _add(self,argument = '', argument2 = '', nargs = None, action = None, default = None, helper = ''):
        self.dc[str(self.name)] = [argument, argument2, nargs, default, action, helper]
        self.name +=1
        return self.dc
    def __repr__(self):
        lines = []
        for k,v in self.dc.items():
            lines.append('{}:{}'.format(k,v))
        return '\n'.join(lines)
