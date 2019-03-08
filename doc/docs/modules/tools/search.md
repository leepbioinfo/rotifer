### Command-line interface package

This package allow you to create a complete command-line interface.

#### Simple usage:

The simplest usage of this package is the following:

```python3
# The first 3 lines the rotifer module path
import os
import sys
sys.insert(0, '/home/kaihami/mymodules')

# Importing the corecli package
import rotifer.core.cli as corecli

# Allow shell autocomplete
import argcomplete

__version__ = 0.01
__authors__ = 'Author1; Author2'

def parse_cli:
   # Create multiple argument parser
   parser = corecli.parser()

   parser.add_argument('--foo') # Add your arguments here

   # Add version control
   parser.add_argument('--version',
        action = 'version',
        version = __version__,
        authors = __authors__,
        program = os.path.basename(__file__),
        description = 'A small description'
                      )

   # Add configdump, configfile and fun
   parser2 = corecli.config().input()

   parser_merged = corecli.merge_parser( parents = [parser, parser2]) # Merge all parsers in one

   argcomplete.autocomplete(parser_merged) # Allow tab autocomplete

   args = corecli.parseargs(parser_merged)

   return args
```

### Actions:

#### Default actions:
Argparse offers a variety of default actions. The most important actions are covered here, for more information please check argparse manual.

#####`append`:
This action stores a list, and append an argument to the list.
```python3
import rotifer.core.cli as corecli

parser = corecli.parser()

parser.add_argument('--mylist',
                    action = 'append')

parser.parse_arg('--mylist apple --mylist banana')

>> Namespace(mylist=['apple', 'banana'])
```

##### **count**
This action count the number of times a keyword occurred. Usefull for verbose control.

```python3
import rotifer.core.cli as corecli

parser = corecli.parser()
parser.add_argument('-v', '--verbose'
                    action = 'count')
```
If the user pass in the command-line -v or --verbose it will count +1. The user can use `-v -v` or `-vv` both will result in verbose = 2.


#### Custom actions:
There are two custom actions, the autoload and autoopen described below.

##### **autoload**
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
    parser.add_argument('accessions',
                        action = corecli.autoload,
                        nargs = '*')

    args = corecli.parseargs(parser)

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
    parser.add_argument('fasta',
                        nargs = '*',
                        action = corecli.autoload,
                        duplicates = False,
                        )

    args = corecli.parseargs(parser)

    return args

if __name__ == '__main__':
    args = parse_cli()
    print(args.fasta)
```


##### `autoopen`
TODO


### Functions:

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
