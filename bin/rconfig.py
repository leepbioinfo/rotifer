#!/usr/bin/env python3
import yaml
import sys
sys.path.insert(0, '/home/kaihami/mymodules')
import argparse
import rotifer.core.cli as corecli
from os.path import expanduser
import os
import sys
import argcomplete
from prettytable import PrettyTable
from datetime import datetime as dt

__version__ = '0.42'
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = argparse.ArgumentParser(description = 'Set configuration files\n\
Usage: rconfig [program] [options]',
                                    formatter_class = argparse.RawTextHelpFormatter)
    subparser = parser.add_subparsers()

    ## pfetch
    pfetch_parser = subparser.add_parser('pfetch')
    pfetch_parser.add_argument('-s', '--show',
                               action = 'store_true')
    pfetch_parser.add_argument('-a','--add',
                               help = 'Add to local database (key:path)')
    pfetch_parser.add_argument('-r','--remove',
                               help = 'Remove from local database (key)')

    ## data
    data_parser = subparser.add_parser('data')

    data_parser.add_argument('-s', '--show',
                               action = 'store_true',
                             help = 'Show data of a project')

    data_parser.add_argument('--show_projects',
                             '-sp',
                             action = 'store_true',
                             help = 'Show projects')

    data_parser.add_argument('-p', '--project',
                             help = 'Add to project')

    data_parser.add_argument('-a','--add',
                               help = 'Add to local database (key:path)')

    data_parser.add_argument('-r','--remove',
                               help = 'Remove from local database (key)')

    data_parser.add_argument('-of','--output_format',
                             help = 'Output format (default: pretty) [pretty/table]',
                             default = 'pretty',
                             dest = 'of')

    #api_key
    api_parser = subparser.add_parser('api_key')
    api_parser.add_argument('-s', '--show',
                               action = 'store_true')
    api_parser.add_argument('-a','--add',
                               help = 'Add to local database (user:api_key)')
    api_parser.add_argument('-r','--remove',
                               help = 'Remove from local database (user)')

    ## acc2operon
    acc_parser= subparser.add_parser('acc2operon')
    acc_parser.add_argument('-s','--show',
                            help = 'Not working',
                               default = 'pfetch.config')

    acc_parser.add_argument('--add',
                            help = 'Not working')
    acc_parser.add_argument('--remove',
                            help = 'Not working')

    create_local_parser = subparser.add_parser('create')
    create_local_parser.add_argument('--new')
    create_local_parser.add_argument('--show_path')
    parser.add_argument('--version',
                        action = 'version',
                        version = corecli.version(program = 'rconfig',
                                                  version = __version__,
                                                  authors = __authors__,
                                                  description = ''))


    create_local_parser.set_defaults(func='create')
    data_parser.set_defaults(func = 'data')
    pfetch_parser.set_defaults(func='pfetch')
    api_parser.set_defaults(func='api_key')

    acc_parser.set_defaults(func='acc2operon')

    argcomplete.autocomplete(parser)

    args = parser.parse_args()

    return args


def showproject():
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config/projects')
    for f in os.listdir(db_local_path):
        print(f.replace('.yaml',''))


def showdata(project, of = 'pretty'):
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config/projects')

    try:
        n = open(os.path.join(db_local_path, project+'.yaml'))
        db_local = yaml.load(n)

    except:
        pass

    try:
        print('## Path: {0}'.format(os.path.join(db_local_path, project+'.yaml')))
        if of == 'pretty':
            table = PrettyTable()
            table.field_names = ['Name', 'Value']
            for k,v in db_local.items():
                table.add_row([k,v])
            print(table)
        else:
            print('Name\tValue')
            for k,v in db_local.items():
                print(k+'\t'+v)


    except:
        print('No project called {0}'.format(project))

def show(db, key_name ,print_all = True):
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config')

    try:
        n = open(os.path.join(db_local_path, db+'.config'))
        db_local = yaml.load(n)

    except:
        pass

    db_global_path = '/home/kaihami/mymodules/rotifer/config/' + db + '.config'
    db_global = yaml.load(open(db_global_path))

    try:
        print('## Local Database')
        print('## Path: {0}'.format(os.path.join(db_local_path, db+'.config')))
        print('Name\tValue')
        for k,v in db_local[key_name].items():
            print(k + '\t' + v)

    except:
        pass
    if print_all:
        try:
            print()
            print('## Global Database')
            print('## Path: {0}'.format(db_global_path))
            print('Name\tValue')
            for k,v in db_global[key_name].items():
                print(k + '\t' + v)
            sys.exit()

        except:
            pass
    sys.exit()


def adddata(s, project, fi):
    k,v = s.split(':')
    dc = {k:v}
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config/projects')

    try:
        db_local = yaml.load(open(os.path.join(db_local_path, project + '.yaml')))
        for k, v in dc.items():
            db_local[k] = v
        with open(os.path.join(db_local_path, project + '.yaml'), 'w') as yaml_file:
            yaml.dump(db_local, yaml_file, default_flow_style=False)

    except:
        try:
            os.makedirs(db_local_path)
        except:
            pass
        for k, v in dc.items():
            db_local = {k:v}
            date_time = {'date': dt.now().strftime('%D')}

        with open(os.path.join(db_local_path, project + '.yaml'), 'w') as yaml_file:
            yaml.dump(db_local, yaml_file, default_flow_style=False)
    print('## Added to {2} [{0}: {1}]'.format(k,v, project))
    showdata(project)
    sys.exit(0)

def add2db(s,key_name , fi):
    k,v = s.split(':')
    dc = {k:v}
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config')
    try:
        db_local = yaml.load(open(os.path.join(db_local_path, fi+'.config')))
        for k, v in dc.items():
            db_local[key_name][k] = v
        with open(os.path.join(db_local_path, fi+'.config'), 'w') as yaml_file:
            yaml.dump(db_local, yaml_file, default_flow_style=False)
    except:
        try:
            os.makedirs(db_local_path)
        except:
            pass

        for k, v in dc.items():
            db_local = {key_name:{k:v}}
        with open(os.path.join(db_local_path, fi+'.config'), 'w') as yaml_file:
            yaml.dump(db_local, yaml_file, default_flow_style=False)
    print('## Added to {2} {0}: {1}'.format(k,v, key_name))
    print()
    show(fi, key_name, print_all = False)
    sys.exit(0)

def removedata(s, project, fi):
    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config/projects')
    try:
        dc = {}
        db_local = yaml.load(open(os.path.join(db_local_path, project+'.yaml')))
        try:
            for k,v in db_local.items():
                print(k, s)
                if k != s:
                    dc[k] = v
            with open(os.path.join(db_local_path, project+'.yaml'), 'w') as yaml_file:
                yaml.dump(dc, yaml_file, default_flow_style=False)
            print('## Removed from {1} {0}'.format(s, project))
            showdata(project)

        except:
            print('No {0} in DB'.format(s))
        sys.exit(0)
    except:
        pass

def remove2db(s,key_name,fi):

    home = expanduser("~")
    db_local_path = os.path.join(home, '.rotifer/config')
    try:
        dc = {key_name:{}}
        db_local = yaml.load(open(os.path.join(db_local_path, fi+'.config')))
        if s in db_local[key_name].keys():
            for k,v in db_local[key_name].items():
                if k != s:
                    dc[key_name][k] = v
            with open(os.path.join(db_local_path, fi+'.config'), 'w') as yaml_file:
                yaml.dump(dc, yaml_file, default_flow_style=False)
            print('## Removed from {1} {0}'.format(s, key_name))
            print()
            show(fi,key_name, print_all = False)
        else:
            print('No {0} in DB'.format(s))
        sys.exit(0)
    except:
        pass

if __name__ == '__main__':

    args = parse_cli()
    try:
        if args.func == 'pfetch':
            if args.show:
                show(args.func, 'database')
                sys.exit()
            elif args.add:
                add2db(args.add,'database', args.func)
                sys.exit()
            elif args.remove:
                remove2db(args.remove, 'database', args.func)
                sys.exit()
            else:
                print('use rconfig pfetch -h/--help for help')

        elif args.func == 'api_key':
            if args.show:
                show(args.func, 'api_key')
                sys.exit()
            elif args.add:
                add2db(args.add, 'api_key', args.func)
                sys.exit()
            elif args.remove:
                remove2db(args.remove, 'api_key', args.func)
                sys.exit()
            else:
                print('use rconfig api_key -h/--help for help')

        elif args.func == 'data':

            if args.show_projects:
                showproject()
                sys.exit()

            elif args.show:
                showdata(args.project, of = args.of)
            elif args.add:
                adddata(args.add, args.project, args.func)
            elif args.remove:
                removedata(args.remove, args.project, args.func)
            else:
                print('use rconfig data -h/--help for help')
                sys.exit()
        if args.func == 'acc2operon':
            pass
    except:
        pass
        # print('use -h/--help for help')
