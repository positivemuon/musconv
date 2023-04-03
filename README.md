# musConv

musConv is a python package for generating near cubic supercells and for checking the convergence of a supercell siz with respect to atomic forces induced by an interstitial impurity at a Voronoi site.

|     | |
|-----|----------------------------------------------------------------------------|
|Latest release| [![PyPI version](https://badge.fury.io/py/musConv.svg)](https://badge.fury.io/py/musConv) [![License](https://img.shields.io/github/license/positivemuon/musConv.svg)](https://pypi.org/project/musConv/) |

## Installation
1.) Install from [pypi](https://pypi.org/project/musConv/0.0.1/) as;

```pip install musConv```


2.) Install the repository as:

```
git clone https://github.com/positivemuon/aiida-musConv.git
cd musConv/
python setup.py install
```



## Available Scripts
```
musConv/
    ├── __init__.py
    └── supcgen.py
    └── chkconv.py
```

## i.) supcgen.py 

Generates a nearly cubic supercell (SC) using the pymatgens [CubicSupercellTransformation](https://pymatgen.org/pymatgen.transformations.advanced_transformations.html). Inserts an intersitial atom (default is hydrogen) in the supercell 
at a Voronoi interstitial site. One of it methods initializes the  supercell generation and the other re-initializes generation of a 
larger supercell-size than the former.

To quickly run the script if installed try; 

 ```python examples/run_supcgen.py```

when not installed try:

```python musConv/supcgen.py examples/LiF.cif```


## ii.) chkconv.py 

This script checks if a supercell (SC) size is converged for muon site calculations
using results of unrelaxed atomic forces from a one shot DFT SCF calculation 
or other potential. Structure input is an ase Atom data while forces as array data.

To quickly run the script if installed try;

```python examples/run_chkconv.py```

or

```python examples/run_chkconv2.py  examples/LiF.cif```

when not installed try:

```python musConv/chkconv.py examples/LiF_p1.cif examples/LiF_p1_forces.txt```

