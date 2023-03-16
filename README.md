# musConv


#1-->supcgen.py
Generates a nearly cubic supercell (SC) for convergence checks.
Inserts a muon in the supercell at a Voronoi interstitial site.
One of it methods initializes the supercell generation and the other 
re-initializes generation of a larger supercell-size than the former.

To quickly run the code try:
```python musConv/supcgen.py examples/LiF.cif```


#2-->chkconv.py
Checks if a supercell (SC) size is converged for muon site calculations
using results of atomic forces from a one shot SCF calculation.

To quickly run the code try:
```python musConv/chkconv.py examples/LiF_p1.cif examples/LiF_p1_forces.txt```


#3--> musConv package can be downloaded and installed  using 
```python setup.py install```

#4 To run the the aiida-musConvworkschain, basic knowlwdge of aiida and requisite aiida-core, plugin  installation and setups are required. How to fix inputs for the workchain can be found in the workchain script. The workchain can be imported after successful installation.


#TO DO
iii) Do documentation
