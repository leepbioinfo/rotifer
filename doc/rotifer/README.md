Rotifer's documentation
=======================

Developtment
------------

### Sequence and alignment modules and programs

- Python
  - ~ Implemented
    - raln2color
      - Feature requests
        - User-defined colors for annotated columns
          - Example: raln2color -a table.tsv (color, name, start, end)
      - Bugs / issues
        - Fix zero-based scale output
    - alnposition
    - sequence\_parser
    - rconsensus
    - aln2fig
      - Expected for library: aln.write("file.png")
  - Desired
    - Core alignment modules / libraries
      - alignment viewer for IPython (like raln2color)
