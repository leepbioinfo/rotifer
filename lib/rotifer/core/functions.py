#!/usr/bin/env python3

from datetime import datetime as dt
import os
import sys
from os.path import expanduser
import yaml
import warnings
import signal
import inspect
from rotifer.core import GlobalConfig
warnings.filterwarnings("ignore")

"""
Rotifer core functions
"""

__version__ = 0.15

def vmsg(message = ''):
    """
    Output a message to stderr.

    Three options available:
    1-) Input a string
    2-) Input a list with one element
    3-) Input a list with n elements

    If a list with more than one element is passed,
    the first element will be printed with timelog the other lines will be printed in a new line.
    -----------
    INPUT:
    message: A message in the format, str or list.
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
    import rotifer.core.functions as cf

    if verbose:
        cd.msg('Small test')
    >> ## [12/03/18 09:50:54] Small test
    """

    now = dt.now().strftime('[%D %H:%M:%S]')
    if isinstance(message, list):
        if len(message) == 0:
            message = ''.join(message)
        else:
            m = message[0] +'\n'
            m2 = ''
            for e in message[1:]:
                m2 += '### ' + e + '\n'
            message = m+m2.rstrip('\n')
    sys.stderr.write('## {0} {1}\n'.format(now, message))
    sys.stderr.flush()

class rdoc:
    '''
    Input a file with rotifer corecli
    Ouput a markdown file
    The first parameter f is the script path
    the md object is created containing a raw markdown
    Later it can be printed to stdout using writer() function.
    ---------
    EXAMPLE:
    import rotifer.core.functions as cf
    md = cf.rdoc(os.path.basename(__file__))
    # Raw markdown
    print(md)

    # Pretty markdown
    md.writer()
    '''

    def __init__(self, f):
        self.f = f
        self.md = self._markdown()

    def _markdown(self):
        __version__ = 0
        __authors__ = ''
        self.f =open(self.f).read().splitlines()
        s2 = []
        for line in self.f:
            if line.startswith('from') or line.startswith('import'):
                s2.append(line)

        docs = []
        for x in range(len(self.f)):
            if self.f[x].startswith('__doc__'):
                for y in range(x+1, len(self.f)):
                    if "'''" in self.f[y]:
                        break
                    else:
                        docs.append(self.f[y])

        dc = {}
        for ele in docs:
            if ele.isupper():
                k = ele
                dc[k] = ''
            else:
                dc[k] += ele+'\n'

        s = []
        for x in range(len(self.f)):
            if 'def parse_cli()' in self.f[x]:
                while True:
                    if self.f[x].lstrip().startswith('return') or self.f[x].lstrip().startswith('args'):
                        break
                    else:
                        x+=1
                        s.append(self.f[x].strip())

        fi = exec('\n'.join(s2), locals(), globals())
        st = '\n'.join(s[:-1])
        exec(st, locals(), globals())
        ls_info = []

        for e in parser_merged._actions:
            long_alias = []
            short_alias = []
            default = []
            arg_type = []

            if e.default != '==SUPPRESS==' and \
                '--fun' not in e.option_strings and \
                '--configdump' not in e.option_strings and \
                '--configfile' not in e.option_strings:
                if not e.option_strings: # empty alias
                    long_alias.append(e.dest)
                    short_alias.append('')
                    default = e.default if not isinstance( e.default, type(None)) else 'None'
                    arg_type = e.type if not isinstance(e.type, type(None)) else 'None'
                    helper = e.help if not isinstance(e.type, type(None)) else ''
                    ls_info.append([long_alias,
                                    short_alias,
                                    default,
                                    arg_type,
                                    helper])

                else:
                    longest_opt = max([len(x) for x in e.option_strings])
                    for option in e.option_strings:
                        if len(option) == longest_opt or len(option)+1 == longest_opt:
                            long_alias.append(option)
                        else:
                            if option:
                                s_option = option
                            else:
                                s_option = ''
                            short_alias.append(s_option)
                    default = e.default if not isinstance( e.default, type(None)) else 'None'
                    arg_type = e.type if not isinstance(e.type, type(None)) else 'None'
                    helper = e.help if not isinstance(e.help, type(None)) else ''
                    ls_info.append([long_alias,
                                    short_alias,
                                    default,
                                    arg_type,
                                    helper])

        toprint = '## Rotifer\n\n'
        head_doc = ['NAME', 'SYNOPSIS', 'DESCRIPTION', 'AUTHOR']
        for descriptor in head_doc:
            if descriptor in dc.keys():
                toprint += '## {0}\n'.format(descriptor)
                toprint += dc[descriptor] +'\n'
            else:
                toprint += '## {0}\n'.format(descriptor)

        toprint += '''
## OPTIONS
'''
        for e in ls_info:
            aliases = ' '.join(e[0]) +', ' + ' '.join(e[1])
            aliases_msg = '**'+ aliases.rstrip(', ')+'**'+'\n'
            default_msg = 'default: '+ str(e[2]) +'\n'
            type_msg = e[3]
            helper = e[4]
            toprint += '\n'.join([aliases_msg, helper, '',default_msg, ''])


        toprint += '''
## Program options summary


|**Long name**|**Aliases**|**Type**|
|-----:|----:|:-----:|
'''

        for e in ls_info:
            long_name = ' '.join(e[0])
            short_name = ' '.join(e[1])
            # type_name = e[3]
            type_name = str(e[3].__class__.__name__) if str(e[3].__class__.__name__) != 'type' else 'None'
            toprint += '|{0}|{1}|{2}|\n'.format(long_name,short_name,type_name)

        toprint = toprint.rstrip('\n')
        return toprint

    def writer(self, output= 'markdown', theme = '641.4291'):
        '''
        '''
        output = self._checker(output)
        self.theme = theme
        switcher = {'md': self._mdwriter,
                    'raw': self._rawwriter,
                    'man': self._manwriter}
        switcher[output]()

    def _checker(self, output):
        if output in ['md', 'markdown']:
            return 'md'
        if output in ['raw']:
            return 'raw'
        if output in ['man']:
            return 'man'
    def _mdwriter(self):
        mdv.term_columns = 60
        print(mdv.main(self.md, theme = self.theme))

    def _rawwriter(self):
        print(self.md)
    def _manwriter(self):
        pass

class dotdict(dict):
    """A class implementing dot annotation access to dictionary attributes
    Instead of using dc['key'] dictionary values can be accessed using dc.key.
    -----------
    PARAMETERS:
    input a python dictionary
    -----------
    Example:
    import rotifer.core.functions as cf

    dc = dotdict({'key1': 'value1',
                  'key2': 'value2'})

    dc.key1
    dc.key2
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def _flatten(*args):
    '''
    Flat a list of arguments
    Up to two levels
    list of list [[]]
    ---------
    PARAMETERS:
    args: list, or singleton
        a list of input
    ---------
    RETURNS:
    A flatten list
    --------
    EXAMPLE:
    args = a
    >> [a]
    args = [a,b]
    >> [a,b]
    args = [[a],b]
    >> [a,b]
    args = [[a,b], [c]]
    >> [a,b,c]
    '''
    ls = []
    if isinstance(args, (list, tuple)):
        for ele1 in args:
            if isinstance(ele1, (list, tuple)):
                for ele2 in ele1:
                    if isinstance(ele2, (list, tuple)):
                        for ele3 in ele2:
                            if isinstance(ele3, (list, tuple)):
                                for ele4 in ele3:
                                    ls.append(ele4)
                            else:
                                ls.append(ele3)
                    else:
                        ls.append(ele2)
            else:
                ls.append(ele1)
    else:
        ls.append(args)
    return ls

def _fasta_ls(*args):
    '''
    Read fasta file or fasta sequence and output a list, each element is a line in the list.
    -----------
    PARAMETERS:
    args: a list containing a fasta sequence or file
    -----------
    RETURNS:
    List: each element is a line of the fasta file
    -----------
    EXAMPLE:

    '''
    fasta_res = []
    for arg in args:
        if isinstance(arg, (list, tuple)):
            fasta_res.extend(list(arg))
        else:
            try:
                fasta_res += open(arg).read().splitlines()
            except:
                fasta_res += arg.split('\n')
    return fasta_res

def not_kwargs(dict_args, key, value):
    """ Small function to check if an element is not in a dictionary.
    It is useful to set a default value to functions that exceeds the
    number of 5 arguments
    """
    if key in dict_args.keys():
        return dict_args[key]
    return value

def loadAPI(username = ''):
    '''
    Load NCBI API key
    '''
    home = expanduser('~')
    db_local_path = None
    for configdir in ["etc","config"]:
        for configext in ["yml","yaml","config"]:
            db_local_path = os.path.join(home, '.rotifer', configdir, f'api_key.{configext}')
            if os.path.exists(db_local_path):
                break
            else:
                db_local_path = None
        if db_local_path:
            break
    db_local = yaml.load(open(db_local_path))

    api_key = ''
    try:
        api_dict = db_local['api_key']
        if username == '':
            api_key = api_dict[list(api_dict.keys())[0]]

        if username in api_dict.keys():
            api_key = api_dict[username]

    except:
        pass
    return api_key

def findDataFiles(load=[], user_path=os.path.join(GlobalConfig['user'],'share'), system_path=os.path.join(GlobalConfig['base'],'share/rotifer/')):
    '''
    This routine locates data files under Rotifer's standard locations.

    User-supplied qualified paths, such as load='./config.yml', are
    accepted but this method is mostly useful when loading files from
    the following standard locations:

      1) The user folder ~/.rotifer/share
      2) <rotifer installation path>/share/rotifer
      3) The path set by the environment variable ROTIFER_DATA

    If the load parameter starts with the prefix ':' and subdirectory
    names are separated by dots ('.'), files will be loaded from the
    standard locations.
    
    Example:
    
      load=":taxonomy.taxonomy.txt"
      
    means
    
      ~/rotifer/share/data/taxonomy/taxonomy.txt

    if the user location above exists, or

      <path to rotifer>/share/rotifer/taxonomy/taxonomy.txt
    
    if only the file from the rotifer installation directory exists.

    Notice that files from the user standard location take precedence
    over files under rotifer's installation path.
    '''

    # Make sure input is a list
    if not isinstance(load,list):
        load = [ load ]

    files = []
    for f in load:
        # Standard files
        if f[0] != ":" and os.path.exists(f):
                files.append(f)

        # Process path
        g = []
        f = f.replace(':','')
        if os.path.sep in f:
            f = f.split(os.path.sep)
        else:
            f = f.split('.')
            if len(f) > 1:
                g = f
                g[-2] = g[-2] + '.' + g[-1]
                del(g[-1])

        # User path
        if os.path.exists(os.path.join(user_path, *f)):
            files.append(os.path.join(user_path, *f))
        elif len(g) > 0 and os.path.exists(os.path.join(user_path, *g)):
            files.append(os.path.join(user_path, *g))

        # System path
        if os.path.exists(os.path.join(system_path, *f)):
            files.append(os.path.join(system_path, *f))
        elif len(g) > 0 and os.path.exists(os.path.join(system_path, *g)):
            files.append(os.path.join(system_path, *g))

        # Databases
        if 'ROTIFER_DATA' in os.environ and os.path.exists(os.path.join(os.environ['ROTIFER_DATA'], *f)):
            files.append(os.path.join(os.environ['ROTIFER_DATA'], *f))
        elif len(g) > 0 and 'ROTIFER_DATA' in os.environ and os.path.exists(os.path.join(os.environ['ROTIFER_DATA'], *g)):
            files.append(os.path.join(os.environ['ROTIFER_DATA'], *g))

    # Return
    return files

def loadConfig(load, user_path=GlobalConfig['userConfig'], system_path=GlobalConfig['baseConfig']):
    '''
    This routine loads configuration parameters from YAML file(s).

    If the load parameter starts with ':' and subdirectory names
    are separated by dots ('.'), the YAML file is searched in the
    user_path directory or, secondly, in the system_path directory.

    Example:
    
      conf = loadConfig(":db.pfam")

    means

      <user_path>/db/pfam.yml

    if <user_path>/db/pfam.yml exists.

    Therefore, if user_path is "~/.rotifer/etc", this
    example translates to

      conf = "~/rotifer/etc/db/pfam.yaml"

    If the location above doesn't exist, the path

      <system_path>/db/pfam.yml
    
    is used, i.e. the user path always takes precedence over
    Rotifer's installation path so that user customizations
    always apply first.
    '''

    if load[0] != ":":
        if os.path.exists(load):
            try:
                return yaml.load(open(load))
            except:
                print(f'Error while parsing file {load}: {sys.exc_info[1]}', file=sys.stderr)
                return {}
        else:
            print(f'File {load} not found.', file=sys.stderr)
            return {}

    # Loading configuration from 
    expand_load = load.replace(':','')

    # User path
    config = {}
    if os.path.exists(user_path):
        try:
            config = yaml_search(expand_load, user_path)
        except:
            print(f'Error while parsing file {os.path.join(user_path,expand_load)}: {sys.exc_info[1]}', file=sys.stderr)

    # System path
    if not config and os.path.exists(system_path):
        try:
            config = yaml_search(expand_load, system_path)
        except:
            print(f'Error while parsing file {os.path.join(system_path,expand_load)}: {sys.exc_info[1]}', file=sys.stderr)

    # Return
    return config

def yaml_search(string, path):
    exts = ['.yml','.yaml','.config']

    data = {}
    ls = string.split('.')
    while ls:

        # Check if any file or directory matches
        if os.path.isdir(path):
            contents = os.listdir(path)
            targets = [ ls[0] + x for x in exts ]
            targets.append(ls[0])
            for target in targets: # Search order: exts, ls[0] (directory)
                if target in contents:
                    path = os.path.join(path, target)
                    break

        # Load matching file that matches an extension
        if os.path.isfile(path) and [x for x in exts if path.endswith(x)]:
            if not data:
                data = yaml.load(open(path))
            # Search for matching elements in data
            if isinstance(data, dict):
                if ls[0] in data.keys():
                    data = data[ls[0]]
            elif isinstance(data,list):
                try:
                    data = data[int(ls[0])]
                except:
                    pass

        # Move stack
        del ls[0]
    #while ls:

    # Return
    return data

def loadClasses(load, user_path=os.path.join(GlobalConfig['user'],'lib'), system_path=os.path.join(GlobalConfig['base'],'lib','rotifer')):
    '''
    This function list all classes and methods
    The main problem is the import
    '''

    classes = {}

    search_file = ''
    yml = ''

    old = os.getcwd()

    if load.startswith(':'):
        yml = loadConfig(load, user_path = user_path,
                         system_path = system_path)

    if yml:
        for path in yml:
            if not search_file:
                os.chdir(path)
                sys.path.insert(0, path)
                for e in os.listdir(path):
                    if e.endswith('.py'):
                        modname = e[:-3]
                        sub = __import__(modname)

                        methods_and_classes = [x for x in dir(sub) if not x.startswith('__')]
                        for method in methods_and_classes:
                            if inspect.isclass(getattr(sub, method)):
                                classes[method] = getattr(sub, method)
                            elif inspect.isfunction(getattr(sub, method)):
                                classes[method] = getattr(sub,method)
                            else:
                                pass

    else:
        if os.path.sep in load:
            load = load.replace(os.path.sep, '.')

        ls = load.split('.')
        if ls[0].startswith('~'):
            path = expanduser('~')

        else:
            path = os.path.join(os.path.sep, ls[0])

        ls = ls[1:]

        while ls:
            if os.path.isdir( path):
                for ele in os.listdir(path):
                    if ele.startswith(ls[0]):
                        if os.path.isdir(os.path.join(path,ele)):
                            path = os.path.join(path, ele)
                        else:
                            pass
            elif ls[0]+'.py' in os.listdir(path):
                search_file = (ls[0])
                break

            del ls[0]

        if not search_file:
            os.chdir(path)
            sys.path.insert(0, path)
            for e in os.listdir(path):
                if e.endswith('.py'):
                    modname = e[:-3]
                    sub = __import__(modname)

                    methods_and_classes = [x for x in dir(sub) if not x.startswith('__')]
                    for method in methods_and_classes:
                        if inspect.isclass(getattr(sub, method)):
                            classes[method] = getattr(sub, method)
                        elif inspect.isfunction(getattr(sub, method)):
                            classes[method] = getattr(sub,method)
                        else:
                            pass

        if search_file:
            while ls:
                if ls[0]+'.py' in os.listdir(path):
                    load_class = ls[0]
                    os.chdir(path)
                    sys.path.insert(0, path)
                    sub = __import__(load_class)
                else:
                    try:
                        classes[ls[0]] = getattr(sub, ls[0])
                    except:
                        pass

                del ls[0]

            if not classes:
                methods_and_classes = [x for x in dir(sub) if not x.startswith('__')]
                for method in methods_and_classes:
                    if inspect.isclass(getattr(sub, method)):
                        classes[method] = getattr(sub, method)
                    elif inspect.isfunction(getattr(sub, method)):
                        classes[method] = getattr(sub,method)
                    else:
                        pass

        os.chdir(old)

    return classes

class wait2finish:
    kill_now = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_nicely)
        signal.signal(signal.SIGTERM, self.exit_nicely)

    def exit_nicely(self,signum, frame):
        self.kill_now = True



def openread(largs):
    '''
    Open, read, and close a file. Remove duplicates
    '''

    # Process sys.argv
    myset=set()
    for iohandle in largs:
        try:
            stdi = iohandle.read().splitlines()
        except:
            stdi = ''
        if ':' in iohandle:
            try:
                _ = open(loadConfig(iohandle))
                myset = myset.union(_.read().splitlines())
            except:
                sys.stderr.write(f'File not found: {iohandle}\n')
                sys.exit()
        else:
            try:
                iohandle = open(iohandle, mode="r")
                myset=myset.union(iohandle.read().splitlines())
                iohandle.close
            except:
                if stdi:
                    myset=myset.union(stdi)
                elif iohandle == "-":
                    myset=myset.union(sys.stdin.read().splitlines())
                else:
                    myset.add(iohandle)
    _ = [x for x in myset if x != '']
    return(_)

def import_path(path):
    '''
    Import Python code from any input file, whatever its extension.
    '''
    import importlib
    module_name = os.path.basename(path).replace('-', '_')
    spec = importlib.machinery.SourceFileLoader(module_name, path)
    spec = importlib.util.spec_from_loader(module_name, spec)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    return module
