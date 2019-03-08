#!/usr/bin/env python3

import argparse
import sys
import copy as _copy
import yaml
import os
from termcolor import colored
### Add arguments ###

__version__ = 0.01

def parser(description =None):
    return argparse.ArgumentParser(add_help = False, description = description)

class config:
    def __init__(self):
        fun_msg = '''
                         (\_/)
                         (^_^)
                         (")(")
                    '''
        not_fun_msg = '''
                           (\_/)
                           (. .)
                          C(")(")
                    '''

        self.dc = {
            'configfile': {'long': '--configfile',
                          'type': argparse.FileType(mode='r')},
            'configdump': {'long': '--configdump'},

                   'fun': {'long': '--fun',
                           'default': [fun_msg, not_fun_msg],
                           'action': 'store_true'}
                }
    def input(self, exclude = [],
              input_order = ['configfile', 'configdump', 'fun'],
              add_help = False):
        self.parser = argparse.ArgumentParser(add_help = add_help)
        for choice in input_order:
            if choice not in exclude:
                if choice == 'fun':
                    self.parser.add_argument(self.dc[choice]['long'],
                                             action = 'store_true')
                elif choice == 'configfile':
                    self.parser.add_argument(self.dc[choice]['long'],
                                             type = self.dc[choice]['type'])
                else:
                    self.parser.add_argument(self.dc[choice]['long'],
                                             action = 'store_true')
        return self.parser


class add_arguments:
    def __init__(self):
        self.dc = {}
        self.name = 1
    def add(self, argument = '', argument2 = '',
                nargs = None, const = None, default = None,
                arg_type = None, choices = None, helper = None,
                metavar = None, action = None, *args):
        if action:
            self.dc[str(self.name)] = [argument, argument2,
                                       nargs, const, default,
                                       arg_type, choice, helper,
                                       metavar, action]
        else:
            self.dc[str(self.name)] = [argument, argument2,
                                       nargs, const, default,
                                       arg_type, choice, helper,
                                       metavar, action]
        return self.dc
    def __repr__(self):
        lines = []
        for k,v in self.dc.items:
            lines.append('{}:{}'.format(k,v))
        return '\n'.join(lines)

### General Actions ###

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
        _ = [x for x in myset if x != '']
        return(_)

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
        _ = [x for x in myset if x != '']
        return(_)

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

    class dict_options_dev(argparse._AppendAction):
        '''
        An action to use with nargs = '1',
        will return a list of dict, separator is '='
        Usage in command line
        -arg my=1 -arg yours=2
        Result:
        args.arg = [{my : 1}, {yours : 2}]
        '''
        def __init__(self, option_strings, default = {}, *args, **kwargs):
            super(action.dict_options, self).__init__(option_strings = option_strings,
                                                    *args, **kwargs)
        def __call__(self, parser, namespace, values, option_string = None):
            items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, []))
            _source, _target = values[0].split('=')
            dd = {}
            dd[_source] = _target
            items.append(dd)
            setattr(namespace, self.dest, items)
    class dict_options(argparse.Action):
        '''
        An action to use with nargs = '1',
        will return a list of dict, separator is '='
        Usage in command line
        -arg my=1 -arg yours=2
        Result:
        args.arg = [{my : 1}, {yours : 2}]
        '''
        def __call__(self, parser, namespace, values, option_string = None):
            items = getattr(namespace,self.dest) or []
            print (items)
            _source, _target = values[0].split('=')
            dd = {}
            dd[_source] = _target
            items.append(dd)
            setattr(namespace, self.dest, items)

    class extend(argparse.Action):

        def __call__(self, parser, namespace, values, option_string=None):
            items = getattr(namespace, self.dest) or []
            items.extend(values)
            setattr(namespace, self.dest, items)


###


## Fix yaml
def parseargs(add_help = True,parents = [], exclude_from_dump = []):
    '''
    Combine multiple argparse.ArgumentParser
    '''
    if len(parents) == 1:
        args = argparse.ArgumentParser(parents=parents, add_help= add_help).parse_args()

    else:
        args = argparse.ArgumentParser(parents=parents, add_help = add_help).parse_args()
    if args.configfile:
        # Load from config file
        data = yaml.load(args.configfile)
        delattr(args, 'configfile')
        for key, value in data.items():
            if isinstance(value, list):
                _ = getattr(args, key)
                if _ is None:
                    setattr(args, key, value)
                else:
                    setattr(args, key, value)
            else:
                setattr(args, key, value)
    if args.configdump:
        ''' Will perform config dump
            Need to deal with list of list, list of list of dict and other formats...
            Need to deal with dict of dict, and so on.
        '''
        try:
            for ele in exclude_from_dump:
                if getattr(args, ele):
                    delattr(args, ele)
        except:
            pass
        print(vars(args))
        delattr(args, 'configdump')
        delattr(args, 'configfile')
        delattr(args, 'accession')
        try:
            delattr(args, 'fun')
        except:
            pass
        cwd = os.getcwd()
        f = os.path.join(cwd, (sys.argv[0].split('/')[-1].replace('.py', '')+ '.config.yaml'))
        print((sys.argv[0].replace('.py', '')+ '.config.yaml'))
        print(f)
        with open(f, 'w') as fi:
            for arg in vars(args):
                if getattr(args,arg) is not None:
                    if isinstance(getattr(args,arg), str): # check if string
                        # Write first bloc
                        s = getattr(args,arg)
                        fi.write(str(arg)+': |\n')
                        fi.write(' '*4)
                        fi.write(getattr(args, arg))
                        fi.write('\n')
                    if isinstance(getattr(args, arg), list): # check if list
                        fi.write(str(arg)+':\n')
                        for ele in getattr(args, arg):
                            fi.write(' '*2 + '- ')
                            if isinstance(ele, dict):
                                for k, v in ele.items():
                                    fi.write(k+': '+str(v))
                                    fi.write('\n')
                            else:
                                fi.write(str(ele))
                                fi.write('\n')
                        fi.write('\n')
                    if isinstance(getattr(args, arg), dict): # check if dict
                        fi.write(str(arg)+':\n')
                        for k,v in getattr(args,arg).items():
                            fi.write(' '*3)
                            fi.write(str(k)+': '+str(v))
                            fi.write('\n')
                        fi.write('\n')
                    if isinstance(getattr(args, arg), bool): # Boolean
                        fi.write(str(arg)+'\n')
                        fi.write(' '*4)
                        fi.write(str(getattr(args, arg)).lower() )
                        fi.write('\n')
                    if isinstance(getattr(args,arg), int):
                        print('HAAHA')
                        print(arg)

        quit()
    return args

## Ok
def version(program = '',
            description = '',
            version = '',
            authors = '',
            program_attrs = {'color': 'red',
                             'attrs': ['bold']
            },
            description_attrs = {'color': 'green',
                             'attrs': []
            },
            version_attrs = {'color': 'cyan',
                             'attrs': []
            },
            authors_attrs = {'color': 'cyan',
                             'attrs': []
            }):
    '''
    Pre-formatted version style
    ------
    Parameters:
    program:     PROGRAM NAME       [String]
    version:     VERSION NUMBER     [String]
    authors:     AUTHORS NAME       [String]
    description: SHORT DESCRITPTION [String]

    Format Parameters:
    program_attrs:     FORMAT STYLE [Dictionary]
    description_attrs: FORMAT STYLE [Dictionary]
    version_attrs:     FORMAT STYLE [Dictionary]
    authors_attrs:     FORMAT STYLE [Dictionary]

    Format Syle:
    color: ['grey', 'red', 'green', yellow', 'blue', 'magenta', 'cyan', 'white']
    attrs: ['bold', 'dark', 'underline', 'blink', 'reverse', 'concealed']
    ------
    Usage:
    > parser.add_argument(--version,
                          action = 'version',
                          version = corecli.version(program = 'my program',
                                                    version = '0.002'
                                                    authors = 'Me',
                                                    program_attrs = {color: 'cyan'}) # Change color to red
                         )

    > parser.add_argument(--version,
                          action = 'version',
                          version = corecli.version(program = 'my program',
                                                    version = '0.002'
                                                    authors = 'Me',
                                                    program_attrs = {color: 'cyan',
                                                                    attrs: ['bold', 'underline']}) # Change color to red, bold and underline
                         )
    '''
    dcs = [program_attrs, description_attrs, version_attrs, authors_attrs]
    ls_keys = ['color', 'attrs']
    for dc in dcs:
        if ls_keys != list(dc.keys()):
            missing = [x for x in ls_keys if x not in list(dc.keys())]
            for e in missing:
                if e == 'color':
                    if dc == program_attrs:
                        dc[e] = 'red'
                    if dc == description_attrs:
                        dc[e] = 'green'
                    else:
                        dc[e] == 'cyan'
                if e == 'attrs':
                    if dc == program_attrs:
                        dc[e] = ['bold']
                    else:
                        dc[e] = []
    s = '''
{p}
{d}
Version: {v}
Authors: {a}
'''.format(space = '#'*20,
           p = colored(program +':', program_attrs['color'], attrs = program_attrs['attrs']),
           d = colored(description, description_attrs['color'], attrs = description_attrs['attrs']),
           v = colored(version, version_attrs['color'], attrs = version_attrs['attrs']),
           a = colored(authors, authors_attrs['color'], attrs = authors_attrs['attrs']),
           )
    return s





def _ensure_value(namespace, name, value):
    if getattr(namespace, name, None) is None:
        setattr(namespace, name, value)
    return getattr(namespace, name)


