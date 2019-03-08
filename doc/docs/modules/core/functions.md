### Core functions package

This package contains the core rotifer functions including:

- [dotdict](#dc)
- [_flatten](#flatten)
- [_fasta_ls](#fastals)
- [loadConfig](#loadconfig)
- [loadClasses](#loadclasses)
- [wait2finish](#wait2finish)
- [optimize_df](#optimizedf)
- [openread](#openread)

### <a name = "dc"> 1. dotdict </a>

This is a silly class to access dictionary attributes using dot annotation.

Common dictionary usage in python uses, dict['KEY'] to access the value. Using this class you can access the attribute using dot (.), as dict.KEY.

It may be usefull to save some typing or sometimes is very handy to use this type of annotation.

Example:
```python
import rotifer.core.functions as rcf

dot = rcf.dotdict({'apple': 0,
                   'banana':10,
                   'orange':20})

print(dot.apple)
>> 0
```

Modifying values
```python
import rotifer.core.functions as rcf

dot = rcf.dotdict({'apple': 0,
            'banana':10,
            'orange':20})
dot.apple = 3
print(dot.apple)

>> 3
```
Passing a dict

```python
import rotifer.core.functions as rcf

common_dictionary = {'apple': 0,
                     'banana':10,
                     'orange':20}

dot = rcf.dotdict(common_dictionary)

print(dot.apple)
>> 0
```

### <a name = "flatten" >2. _flatten </a>

This function allow extend and flatten multiple lists and/or strings, up to two levels.

Below is an example using strings as arguments for this function.

```python3
import rotifer.core.functions as rcf

args = 'apple'
flat = rcf._flatten
flat(args)
>> ['apple']

args2 = 'banana'

flat(args,args2)
>> ['apple', 'banana']

flat(args,args2, 'pineapple')
>> ['apple', 'banana', 'pineapple']
```

A more usefull usage is when mixing strings and lists/tuples.

```python3
import rotifer.core.functions as rcf

flat = rcf._flatten

# String
args = 'apple'

flat(args)
>> ['apple']

# A list
args2 = ['banana', 'melon']

flat(args, args2)
>> ['apple', 'banana', 'melon']

# Tuple
args3 = ('pineapple',)
type(args3)
>> tuple

flat(args, args2, args3)
>> ['apple', 'banana', 'melon', 'pineapple']
```

In the above example, notice that if you input a tuple they will be transformed to a flat python list.

This function allows to flatten a `list of list` and/or a `list of list of lists`.

```python3
import rotifer.core.functions as rcf

flat = rcf._flatten

# list of lists
args = [['apple', 'banana']]

flat(args)
>> ['apple', 'banana']

# list of lists
args2 = [['orange', 'pineapple'], 'melon']

# Flatten two list of lists
flat(args, args2)
>> ['apple', 'banana', 'orange', 'pineapple', 'melon']

# The input could be a list of list of lists
args3 = [[('cherry', 'raspberry')]]

flat(args3)
>> ['cherry', 'raspberry']
```

### <a name = "fastals" >3. _fasta_ls</a>

This function accepts a fasta file, fasta in the string format, and/or a fasta in a list format.

```python3
import rotifer.core.functions as rcf

fasta_ls = rcf._fasta_ls
fasta1 = '''>Seq1
AATCGCATTA
>Seq2
TTACCGGGAAA
'''

fasta2 = ['>Seq3', 'ACGTCGA', '>Seq4', 'AAATTTCA']

fasta3 = '/path/to/fasta.file'

fasta_ls(fasta1, fasta2, fasta3)
```

### <a name = "loadconfig" >4. loadConfig</a>
This function allow to load a configuration file in the yaml format.

Important parameters includes:

<center>

|Parameter|description|
|:---:||
|load|A file to be loaded|
|user_path||
|system_path||

</center>


### <a name = "loadclasses" >5. loadClasses</a>
This function allows to load a custom classes from any python file.

<center>

|Parameter|description|
|:-------:|-----------|
|load     | A file to be loaded|
|user_path||
|system_path||

</center>

This is usefull for generic switchers.

Example:

In a file and/or a folder containing multiple files. Each file contains a class or classes to be loaded.

You can create a file named `classes.py`

`classes.py`
```python3
# A file containing multiple classes

class person:
    def __init__(self, name, legs):
        self.name = name
        self.legs = legs
    def writer(self):
        print(f'Hey I'm a person and my name is: {self.name} and I have {self.legs} legs')

class animal:
    def __init__(self, name, legs):
        self.name = name
        self.legs = legs
    def writer(self):
        print(f'Hey I'm an animal and my name is: {self.name} and I have {self.legs} legs')
```

In another file contains the a generic class. This class (named `loader` in this example) has a switch case which will load all classes present in the `classes.py` file, in this case `animal` and `person`.

The loader class accepts an argument, here called source, as well as the two mandatory paramenters for the animal and person classes (name and legs).
The loadClasses function loads all classes (or a specific class if specified) returning a dictionary, the key is the class name and the value in the class itself. One example of how to use loadClasses as a generic loader is showed below, in the self.switcher (creating the dictionary) the class loader behaviour will depends on the source type (animal or person), and the name and legs parameters are passed in the self.target line.

```python3
import rotifer.core.functions as rcf

class loader:
    def __init__(self, source, name, legs):
        self.switcher = rcf.loadClasses('./classes.py') # path to file
        self.target = switcher[source](name, legs)

    def writer(self):
        return self.target.writer()

    def show_sources(self):
        return self.switcher.keys()

test_person = loader(source = 'person', name = 'John', legs = 2)
test_person.writer()
>> Hey I'm a person and my name is: John and I have 2 legs

test_animal = loader(source = 'animal', name = 'Fido', legs = 4)
test_animal.writer()
>> Hey I'm an animal and my name is: Fido and I have 4 legs

print('\n'.join(test_animal.show_sources()))
>> animal
>> person
```

loadClasses can accept a yaml file as well. The yaml file must contain the path where the classes are. This yaml must be stored inside the rotifer module path inside the config folder and/or the user folder (.rotifer/config).

The loadClasses will invoke loadConfig if the argument passed begins with ':' (see [loadConfig](#loadconfig))

For example
```python3
loadClasses(':classfolder.classfile')
```

This yaml file is located inside /rotifer/config/classfoder/ and named as classfile.yaml
```yaml
- /path/to/class
- /another/path/to/class
```

In this example, loadClasses will load all classes located in the two paths insite the yaml file.

### <a name = "wait2finish" >6. loadClasses</a>

### <a name = "optimizedf" >7. optimize_df</a>

### <a name = "openread" >8. openread</a>


