ROTIFER
=======

Rapid Open-source Tools and Infrastructure For data Exploration and Research
----------------------------------------------------------------------------

ROTIFER is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.


Installing ROTIFER package with conda enviroments:

1. Clone this git reposytory.
2. Create a rotifer's conda enviroment:
  ```Bash
  conda create --name rotifer biopython pandas numpy termcolor PyYAML tqdm ipython ascii_graph
  ```
3. Create a rotifer profile to the ipython:
  ```
  mkdir ~/.ipython/profile_rotifer
  mkdir ~/.ipython/profile_rotifer/startup
  
  


