#!/usr/bin/env python3

import argparse
import sys
import copy as _copy
import yaml
import os
sys.path.insert(0, '/home/kaihami/mymodules')
# from rotifer.core.loadpath import path as rpath
from termcolor import colored
from collections import defaultdict
from os.path import expanduser
import warnings

#########
### This is a hack we should think carefully about it
#########
import __main__ as mainprog

warnings.filterwarnings("ignore")

__version__ = 0.15
__authors__ = 'Gilberto Kaihami; Robson Souza'

### Add arguments ###
def check_args(dc):
    '''
    Check if a argument is present or not in a config dict
    '''
    ls_args = ['long_arg', 'short_arg', 'dest',
               'nargs', 'const','default',
               'arg_type',
               'choices',
               'helper',
               'metavar', 'action']

    dc_with_args = dotdict({})

    for arg in ls_args:
        if arg in dc.keys():
            dc_with_args[arg] = dc[arg]
        else:
            dc_with_args[arg] = None

    return dc_with_args

class parser:
    def __init__(self, description = '',
                 add_help = True):
        '''
        Create an ArgumentParser
        PARAMETERS:
        description: The program description (String)
        add_help:    Add helper (Boolean)
        -----------
        OUTPUT:
        A parser object
        -----------
        USAGE:
        import rotifer.core.cli as corecli
        parser = corecli.parser()
        parser.add(ARGUMENT)
        args = parser.parse_args()

        return args
        '''

        self.parser = argparse.ArgumentParser(description = description,
                                              formatter_class = argparse.RawTextHelpFormatter,
                                              add_help = add_help)

        self._action_switch = {'action.autoload': action.autoload,
                          'action.autoopen': action.autoopen,
                          'action.fun': action.fun
                         }

        self._type_switch = {'read': argparse.FileType(mode = 'r')
                        }

        self._helper_switch = {'supress':argparse.SUPPRESS
                               }

    def add(self, long_arg = None, short_arg = None, dest = None,
            nargs = None, const = None, default = None,
            arg_type = None, choices = None, helper = None,
            metavar = None, action = 'store', rconfig = '',
            version = '',
            **kw):
        '''
        Add an argument to parser
        ------------
        PARAMETERS:
        long_arg:  long alias
        short_arg: short alias
        dest:      value destination
        nargs:     Number of args
        const:     see argparse const
        default:   Default value
        arg_type:  argument type (no quotes)
        choices:   offer a list of choices
        helper:    Help message
        metavar:   Metavar
        action:    An action
        rconfig:   Load pre-defined configurations
        ------------
        USAGE:
        Add your own argument:
        parser.add(long_arg = '--long',
                   short_arg = '-short',
                   dest = 'myvar'
                   )

        parser.add(':cli.core') # load predefined configuration
        parser.parse_args()
        '''

        # Parse rconfig
        if long_arg is not None:
            if ':' in long_arg:
                rconfig = long_arg
                long_arg = None
            else: pass

        if rconfig:
            home = expanduser("~")
            db_local_path = os.path.join(home, '.rotifer/config')
            f = '/'.join(rconfig.lstrip(':').split('.'))+'.yaml'
            arg2parse = ''
            db_local = ''

            try:
                db_local = yaml.load(open(os.path.join(db_local_path, f)))
            except:
                pass

            try:
                db_global_path = '/home/kaihami/mymodules/rotifer/src'
                db_global = yaml.load(open(os.path.join(db_global_path, f)))
            except:
                pass

            if db_local:
                arg2parse = db_local
            else:
                arg2parse = db_global

            for k in arg2parse.keys():
                parsed_args = check_args(arg2parse[k])

                # Parse custom action
                if parsed_args.action in self._action_switch.keys():
                    parsed_args.action = self._action_switch[parsed_args.action]
                else:
                    pass

                # Parse custom type
                if parsed_args.type in self._type_switch.keys():
                    parsed_args.type = self._type_switch[parsed_args.type]
                else: pass

                # Parse helper
                if parsed_args.helper in self._helper_switch.keys():
                    parsed_args.helper = self._helper_switch[parsed_args.helper]

                if parsed_args.action in ['store_true', 'store_false']:
                    if parsed_args.short_arg is None and parsed_args.long_arg is None:
                        self.parser.add_argument(dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action,
                                                 kw = parsed_args.kw
                                                 )

                    if parsed_args.short_arg is not None:
                        self.parser.add_argument(parsed_args.long_arg,
                                                 parsed_args.short_arg,
                                                 dest = parsed_args.dest,
                                                 help = parsed_args.helper,
                                                 action = parsed_args.action
                                                 )

                    else:
                        self.parser.add_argument(parsed_args.long_arg,
                                                 # parsed_args.short_arg,
                                                 dest = parsed_args.dest,
                                                 help = parsed_args.helper,
                                                 action = parsed_args.action
                                                 )

                elif 'kw' in arg2parse[k]:
                    if parsed_args.short_arg is None and parsed_args.long_arg is None:
                        self.parser.add_argument(dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action,
                                                 kw = parsed_args.kw
                                                 )

                    elif parsed_args.short_arg is not None:
                        self.parser.add_argument(parsed_args.long_arg,
                                                 parsed_args.short_arg,
                                                 dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action,
                                                 kw = parsed_args.kw
                                                 )

                    else:
                        self.parser.add_argument(parsed_args.long_arg,
                                                 # parsed_args.short_arg,
                                                 dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action,
                                                 kw = parsed_args.kw
                                                 )

                else:
                    if parsed_args.short_arg is None and parsed_args.long_arg is None:
                        self.parser.add_argument(dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action
                                                 )

                    elif parsed_args.short_arg is not None:
                        self.parser.add_argument(parsed_args.long_arg,
                                                 parsed_args.short_arg,
                                                 dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action
                                                 )

                    else:
                        self.parser.add_argument(parsed_args.long_arg,
                                                 # parsed_args.short_arg,
                                                 dest = parsed_args.dest,
                                                 nargs =parsed_args.nargs,
                                                 const = parsed_args.const,
                                                 default = parsed_args.default,
                                                 type = parsed_args.arg_type,
                                                 choices = parsed_args.choices,
                                                 help = parsed_args.helper,
                                                 metavar = parsed_args.metavar,
                                                 action = parsed_args.action
                                                 )


        # if User input
        else:
            if action in ['store_true', 'store_false']:
                if short_arg is not None:
                    self.parser.add_argument(long_arg,
                                             short_arg,
                                             dest = dest,
                                             help = helper,
                                             action = action
                                             )

                else:
                    self.parser.add_argument(long_arg,
                                             # short_arg,
                                             dest = dest,
                                             help = helper,
                                             action = action
                                             )

            elif action in ['version']:
                self.parser.add_argument(long_arg,
                                         action = action,
                                         version = version
                                         )
            elif action in ['count']:
                if short_arg is None:
                    self.parser.add_argument(long_arg,
                                            action = action,
                                            help = helper,
                                            dest = dest)
                else:
                    self.parser.add_argument(long_arg,
                                            short_arg,
                                            action = action,
                                            help = helper,
                                            dest = dest)



            elif kw:
                if short_arg is None and long_arg is None:
                    self.parser.add_argument(dest = dest,
                                             nargs =nargs,
                                             const = const,
                                             default = default,
                                             type = arg_type,
                                             choices = choices,
                                             help = helper,
                                             metavar = metavar,
                                             action = action,
                                             kw = kw
                                             )

                elif short_arg is not None:
                    self.parser.add_argument(long_arg,short_arg,
                                             dest = dest,
                                             nargs =nargs,
                                             const = const,
                                             default = default,
                                             type = arg_type,
                                             choices = choices,
                                             help = helper,
                                             metavar = metavar,
                                             action = action,
                                             kw = kw
                                             )
                else:
                    self.parser.add_argument(long_arg,
                                             dest = dest,
                                             nargs =nargs,
                                             const = const,
                                             default = default,
                                             type = arg_type,
                                             choices = choices,
                                             help = helper,
                                             metavar = metavar,
                                             action = action,
                                             kw = kw
                                             )

            else:
                if short_arg is None and long_arg is None:
                    self.parser.add_argument(dest = dest,
                                             nargs =nargs,
                                             const = const,
                                             default = default,
                                             type = arg_type,
                                             choices = choices,
                                             help = helper,
                                             metavar = metavar,
                                             action = action
                                             )
                elif short_arg is not None:
                    self.parser.add_argument(long_arg,short_arg,
                                             dest = dest,
                                             nargs =nargs,
                                             const = const,
                                             default = default,
                                             type = arg_type,
                                             choices = choices,
                                             help = helper,
                                             metavar = metavar,
                                             action = action
                                             )

                else:
                    self.parser.add_argument(long_arg,
                                             dest = dest,
                                             nargs =nargs,
                                             const = const,
                                             default = default,
                                             type = arg_type,
                                             choices = choices,
                                             help = helper,
                                             metavar = metavar,
                                             action = action
                                             )

        return self.parser

    def parse_args(self, exclude_from_dump = [], exclude_core = False):
        '''
        Parse args
        Exclude one argument from configdump output
        IMPORTANT: Need to fix configfile
        -----------
        INPUT:
        exclude_from_dump: A list containing if you want to exclude from the output
        ----------
        OUTPUT:
        A parsed argument
        ---------
        USAGE:
        Simple usage:
        parser = corecli.parser()
        parser.add(':core')
        args = parser.parse_args()
        return args

        Advanced usage:
        parser = corecli.parser()
        parser.add('file')
        parser.add(':core')
        args = parser.parse_args(exclude_from_dump = ['file']) # File arguments not present in the output!
        return args
        '''
        if not exclude_core:
            self.add(':cli.core')
        else:
            pass
        try:
            mainprog.__version__
            mainprog.__authors__
        except:
            print('## ERROR!')
            print('Include __version__ and/or __authors__ in the main program')
            sys.exit()

        self.add('--version',
                 action = 'version',
                 version = version(program = mainprog.__file__,
                             version = mainprog.__version__,
                             authors = mainprog.__authors__,
                             description = self.parser.description))

        if '--configdump' in sys.argv or '--configfile' in sys.argv:

            if '--configdump' in sys.argv:
                for i, a in enumerate(self.parser._actions):
                    if a.nargs == '*':
                        try:
                            del self.parser._actions[i]
                        except:
                            pass

            args = self.parser.parse_args()

            known_args = [x for x in list(args.__dir__()) if not x.startswith('_')]
            if args.configfile:
                '''
                Strategy: Check args,
                '''

                known_args = [x for x in list(args.__dir__()) if not x.startswith('_')]
                data = yaml.load(open(args.configfile))
                delattr(args, 'configfile')
                parents = vars(args).keys()

                all_defaults = {}
                for key in vars(args):
                    all_defaults[key] = self.parser.get_default(key)

                # Add extend option
                extend = False
                for key, value in data.items():
                    if isinstance(value, list): # Check if value in the config file is a list
                        _ = getattr(args, key)

                        all_ls = []
                        all_ls.extend(value)

                        if _:
                            all_ls.extend(_)
                        else:
                            pass
                        if extend:
                            if _ is not None:
                                setattr(args, key, all_ls)
                            else:
                                setattr(args, key, all_ls)
                        else:
                            if _ is not None and _ != all_defaults[key]:
                                setattr(args, key, _)
                            else:
                                setattr(args, key, value)

                    elif isinstance(value, bool): # Set boolean
                        d_default = all_defaults[key]

                        _ = getattr(args, key)
                        if value != d_default:
                            setattr(args, key, value)

                        elif _ != d_default or _ != value:
                            setattr(args, key, _)

                    else:
                        _ = getattr(args, key)
                        d_default = all_defaults[key]
                        if value != d_default:
                            setattr(args, key, value)

                        elif _ != d_default and _ != value:
                            setattr(args, key, _)

                        else :
                            setattr(args, key, d_default)

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
                try:
                    delattr(args, 'configdump')
                except: pass
                try:
                    delattr(args, 'configfile')
                except: pass

                try:
                    delattr(args, 'accession')
                except: pass
                try:
                    delattr(args, 'fun')
                except:
                    pass
                for arg in vars(args):
                    if getattr(args,arg) is not None:
                        k = {arg: getattr(args, arg)}

                        yaml.dump(k, sys.stdout, default_flow_style = False)# default_flow_style = False))
                quit()
            return args
        else:
            args = self.parser.parse_args()
            return args


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
            duplicates = True, remove duplicate lines
            duplicates = False, keep duplicate lines
        '''
        def __init__(self, option_strings, file_type = 'list', kw = {'duplicates':True},
                      *args, **kwargs):
            duplicates = kw['duplicates']
            self._add1 = file_type
            self._add2 = duplicates
            super(action.autoload, self).__init__(option_strings=option_strings,
                                           *args, **kwargs)
        def __call__(self, parser, namespace, values, option_string = None):
            #TODO insert new function for the var file type
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

    class fun(argparse.Action):
        '''
        Fun or not fun message
        Usage:
        if fun:
            sys.stderr.write(fun.fun)
        or
        if fun:
            # Error?
            sys.stderr.write(fun.notfun)
        '''
        def __init__(self, option_strings, *args, **kwargs):
            super(action.fun, self).__init__(option_strings = option_strings, *args, **kwargs)
            self.nargs = 0
        def __call__(self, parser, namespace, values, option_string = None):
            fun_msg = '''
 (\_/)
 (^_^)
 (")(")
'''
            nofun_msg = '''
 (\_/)
 (. .)
C(")(")
'''
            f = dotdict({'fun':fun_msg,
                         'nofun': nofun_msg})
#            f2 = action().nofun
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

class dotdict(dict):
    "dot.notation access to dictionary attributes"
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


