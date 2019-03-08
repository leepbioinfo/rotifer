### Command-line interface package

This package allow you to create a complete command-line interface.

#### 1. Simple usage:

##### 1.1. Creating the parse_cli

The simplest usage of this package is the following:

```python3
import rotifer.core.cli as corecli

__version__ = 0.001
__authors__ = 'YOUR NAME'

def parse_cli():
    parser = corecli.parser(description = 'PROGRAM DESCRIPTION')

    # Add another options here

    args = parser.parse_args()

    return args
```
Firstly you need to import the core CLI module (`import rotifer.core.cli`).



For every rotifer tool the version control is mandatory (`__version__`) and the authors is optional (`__authors__`), although highly recommended.

In the example above, the parse_cli function will return the parsed arguments from the command-line.

As default, this function will add four arguments:

- --help
- --version
- --configfile
- --configdump

Configdump dumps the command-line arguments to a configuration file, while configfile loads this configuration file to the main program.

#### 2. Adding arguments

##### 2.1. Custom arguments

To add new arguments to the parser use `add` method.

```python3
def parse_cli():
    parser = corecli.parser(description = 'PROGRAM DESCRIPTION')

    # Add another options here

    parser.add( long_arg = '--long',
                short_arg = '-short',
                dest = 'myvar',
                nargs = None,
                default = None,
                arg_type = None,
                helper = 'Some help',
                action = "store"
                )

    args = parser.parse_args()

    return args
```

##### 2.2. Arguments from a configuration file

To add new arguments from a configuration file use the following syntax (:cli.[yaml.file.name]). The ":" in the beginning indicates that the file will be loaded from a yaml file (see [core functions](./functions.md/#loadConfig) for a more detailed explanation).

```python3
def parse_cli():
    parser = corecli.parser(description = 'PROGRAM DESCRIPTION')

    # Add another options here

    parser.add(':cli.acc')

    args = parser.parse_args()

    return args
```

The configuration file must be in the .rotifer/config folder (user folder) or in the system folder (/rotifer/module/path/config/cli).

Currently there are only two configurations available at system level, the core configuration and acc configuration.

As mentioned above, the core configuration file is automatically loaded using the corecli module.

This core configuration file contains:

- --help
- --version
- --configfile
- --configdump

The acc configuration file (loaded using, `parser.add(':cli.acc')`) automatically load as list of accessions from the command-line and/or a file. The values stored is a non-redudant python list (removes duplicated elements from the list), and can be accessed using .accession, as showed below.

```python3
def parse_cli():
    parser = corecli.parser(description = 'PROGRAM DESCRIPTION')

    parser.add(':cli.acc')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()
    print(args.accession)
```

### 3. Actions

Rotifer corecli expands the default actions from `argparse` module. It is important to note that the default argparse actions, as well as the custom actions can be used.

#### 3.1. Default actions

Argparse offers a variety of default actions. The most important actions are covered here, for more information please check argparse manual.

##### 3.1.1 `append`:
This action stores a list, and append an argument to the list.

```python3
import rotifer.core.cli as corecli

parser = corecli.parser()

parser.add('--mylist',
           dest = 'mylist',
           action = 'append')

parser.parse_arg('--mylist apple --mylist banana')

>> Namespace(mylist=['apple', 'banana'])
```

##### 3.1.2 `count`

This action count the number of times a keyword occurred. Usefull for verbose control.

```python3
import rotifer.core.cli as corecli

parser = corecli.parser()
parser.add('-v', '--verbose',
           dest = 'verbose',
           action = 'count')
```

If the user pass in the command-line -v or --verbose it will count +1. The user can use `-v -v` or `-vv` both will result in verbose = 2.

#### 3.2. Custom actions

There are two custom actions, the autoload and autoopen as described below.

##### 3.2.1 `autoload`

This action open and load a file in the memory returning a list of strings.

Also, this action can remove (or not) duplicated lines in the file (non-redundant accession for example) and accepts pipe ( '|' ) from the command-line.

To remove duplicated lines (default) set `duplicates = False`.

To keep duplicated lines set `duplicates = True`.

One important behaviour when `duplicates = False` is that the original order is not the same as the input.

Usually this action is used with `nargs = '*'` as showed in the example below.

Example:
Removing duplicated lines. For example, to remove duplicated accessions from a file
```python3
import rotifer.core.cli as corecli

def parse_cli():
    # Create an argument parser
    parser = corecli.parser()

    # Add argument with autoload action

    parser.add(dest = 'accessions',
               action = corecli.autoload,
               nargs = '*')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()
    print(args.accessions)
```

Keep duplicated lines in the file. For example, when reading a fasta file (each line is an element in the list).

```python3
import rotifer.core.cli as corecli

def parse_cli():
    # Create an argument parser
    parser = corecli.parser()

    # Add argument with autoload action
    parser.add(dest = 'fasta',
               nargs = '*',
               action = corecli.autoload,
               duplicates = False,
               )

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()
    print(args.fasta)
```

##### `autoopen`
TODO

### 4. Methods:

#### parser
Add a argument parser

This function substitutes argparse.ArgumentParser().
If you decide to use this function there is no need to import argparse in the program.
The add_help is set to False. This is handy if you use `corecli.merge_parser` (see above).

```python
def parser(description = None):
    return argparse.ArgumentParser(add_help = False, description = description,
                                   formatter_class = argparse.RawTextHelpFormatter)
```

Usage example:
```python
import rotifer.core.cli as corecli

parser = corecli.parser()
parser.add_argument(--foo) # same syntax as argparse
...

```


#### merge_parser

This function merge multiple parsers in one.

For example, with this function is possible to merge a custom parser with the core parser (containing configfile, configdump and fun).

```python3
import rotifer.core.cli as corecli

def parse_cli:
   # Create multiple argument parser
   parser = corecli.parser()

   parser.add_argument('--foo') # Add your arguments here

   # Add configdump, configfile and fun
   parser2 = corecli.config().input()

   # A new parser
   parser3 = corecli.parser()
   parser3.add_argument('--bar')

   parser_merged = corecli.merge_parser( parents = [parser, parser2, parser3]) # Merge all parsers in one

   args = corecli.parseargs(parser_merged)

   return args
```
In this example, three parsers were created: parser, parser2, and parser3. The arguments were merged using merge_parser() function. Usually this function is usefull to merge the rotifer cli with a program specific parser (see [Simple Usage](#simple-usage))

#### parseargs

This function parse the arguments, it is important to use this function if configdump/configfile is present. It will dump/update the values.

You can exclude arguments to be not present in the config file, usually a user input file. This is possible using exclude_from_dump option, the input is a list with argument to exclude.  
Example:
```python3
import rotifer.core.cli as corecli

def parse_cli():
    # Create an argument parser
    parser = corecli.parser()

    # Add argument with autoload action
    parser.add_argument('accessions',
                        action = corecli.autoload,
                        nargs = '*')

    parser2 = corecli.config().input()

    parser_merged = corecli.merge_parser(parents = [parser, parser2])

    args = corecli.parseargs(parser,
                            exclude_from_dump = ['accessions'])

    return args

if __name__ == '__main__':
    args = parse_cli()
    ...
```

#### A fun option
In the corecli.config().input() there is a hidden option called fun. To activate this flag the user should input `myscript --fun`. This will turn `args.fun == True` and you can print to the user two outputs (fun or nofun).

Usually if the programs finished with no errors you can print the fun message, otherwise the nofun message.

A simple case of fun or nofun is the following

```python
import sys
import rotifer.core.cli as corecli

def parse_cli():
    # Create an argument parser
    parser = corecli.parser()

    # Add argument with autoload action
    parser.add_argument('accessions',
                        action = corecli.autoload,
                        nargs = '*')

    parser2 = corecli.config().input()

    parser_merged = corecli.merge_parser(parents = [parser, parser2])

    args = corecli.parseargs(parser,
                            exclude_from_dump = ['accessions'])

    return args

if __name__ == '__main__':
    args = parse_cli()
    try:
        # Long program
        # finished everything
        if args.fun:
            sys.stderr.write(args.fun.fun)

    except ...: # Catch some error
        if args.fun:
            sys.stderr.write(args.fun.nofun)
```
