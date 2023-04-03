# musConv

|     | |
|-----|----------------------------------------------------------------------------|
|Latest release| [![PyPI version](https://badge.fury.io/py/musConv.svg)](https://badge.fury.io/py/musConv) [![License](https://img.shields.io/github/license/positivemuon/musConv.svg)](https://pypi.org/project/musConv/) |

The musConv package is a python package including an [AiiDA](www.aiida.net) workchain for generating a converged supercell for and interstitial defect calculation e.g for muon calculations. 

The main package contains three scripts (classes), two of which are independent and can be used separately and the third an [AiiDA](www.aiida.net) workchain.

At the moment, the musConv package can be downloaded and installed  using:

```python setup.py install```

or with [pypi](https://pypi.org/project/musConv/0.0.1/) using;

```pip install musConv```

The scripts; what they do and usage:


#*******

#1-->supcgen.py script

Generates a nearly cubic supercell (SC) using the pymatgens [CubicSupercellTransformation](https://pymatgen.org/pymatgen.transformations.advanced_transformations.html).
Inserts an intersitial atom (default is hydrogen) in the supercell 
at a Voronoi interstitial site. One of it methods initializes the 
supercell generation and the other re-initializes generation of a 
larger supercell-size than the former.

To quickly run the script if installed try; 

 ```python examples/run_supcgen.py```

when not installed try:

```python musConv/supcgen.py examples/LiF.cif```



#*******

#2-->chkconv.py script

This script checks if a supercell (SC) size is converged for muon site calculations
using results of unrelaxed atomic forces from a one shot DFT SCF calculation 
or other potential. Structure input is an ase Atom data while forces as array data.

To quickly run the script if installed try;

```python examples/run_chkconv.py```

or

```python examples/run_chkconv2.py  examples/LiF.cif```

when not installed try:

```python musConv/chkconv.py examples/LiF_p1.cif examples/LiF_p1_forces.txt```



****
#TO DO

iii) improve documentation
