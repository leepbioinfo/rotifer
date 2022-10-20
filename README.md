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
  conda create --name rotifer biopython pandas numpy termcolor PyYAML tqdm ipython ascii_graph matplotlib ete3
  ```
3. Activating the rotifer enviroment:
  ```Bash
  conda activate rotifer
  ```
4. Setting up the rotifer tools at the conda rotifer enviroment:
  ```
  mkdir -p $CONDA_PREFIX/etc/conda/activate.d
  mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
  touch ./etc/conda/activate.d/env_vars.sh
  touch ./etc/conda/deactivate.d/env_vars.sh
  conda-develop /path_to_git/rotifer/lib
  
  5. Write the follow path to eh activate.d file:
  export PATH="path_to_git/rotifer/bin:$PATH"
  ```
  
  ## Activating debbug function
 ```python
import rotifer
rotifer.logger.setLevel(rotifer.logging.DEBUG)
  ```


