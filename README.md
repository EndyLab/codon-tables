# codon-tables
Code relating to arbitrary DNA decoding and codon table representation, visualization and optimization.

All code is written in Python 3.5. Below is a brief description of the modules and Jupyter Notebooks included with this repo.

## codonUtils
A module used to package general purpose data and methods for codon table manipulation and representation (e.g. list of dNTPs, list of triplet codons, Polar Requirement Scale, block structure of the standard table, method for getting amino acid degeneracy, etc.). Import as follows:
```python
from src.codonUtils import utils
help(utils)
```
## codonTable
A class used to handle representation and visualization of codon table objects. Can be used to represent a given codon table as a pandas DataFrame, a python dict, a networkx graph and corresponding adjacency matrix (representing the probability of mutating between amino acids given a table), and as plots visualizing these types. Import as follows:
```python
from codonTable import codonTable
help(codonTable)
``` 

## codonOptimizer
A class used to handle codon table space search and optimization. Includes definitions for cost functions as well as MCMC style optimization algorithm(s). Import as follows:
```python
from src.codonOptimizer import tableOptimizer
help(tableOptimizer)
```

## Decoding Notebook
A notebook used to demonstrate the codonTable class

## Network Stats Notebook
A notebook used to calculate and visualize statistics of random and optimized codon tables with respect to the standard table.
