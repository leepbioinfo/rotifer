#!/usr/bin/env python3

import os
import inspect
from os.path import expanduser
import sys
import yaml
from rotifer.core.functions import loadConfig

# def getClasses(load = '',
#         system_path = os.path.join(os.path.realpath(os.path.join( os.path.abspath(__file__).split('rotifer')[0]+'rotifer')))):
#     expand_load = load.replace(':', '')


def get_classes(load,
        user_path = os.path.expanduser(os.path.join(*['~', '.rotifer', 'config'])),
        system_path = os.path.join(os.path.realpath(os.path.join( os.path.abspath(__file__), '..', '..')), 'config')):
    '''
    This function list all classes and methods
    The main problem is
    '''
    # load = 'home.kaihami.mymodules.rotifer.test.dyn.dyn'

    classes = {}
    search_file = ''
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
            print(classes)




    else:
        ls = load.split('.')


        if ls[0].startswith('~'):
            path = expanduser('~')
        else:
            path = os.path.join(os.path.sep, ls[0])

        ls = ls[1:]

        while ls:
            if os.path.isdir( path):

                for ele in os.listdir(path):
                    if ls[0] in ele:
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

