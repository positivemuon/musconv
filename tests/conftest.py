import pytest
from pymatgen.core import Structure
import numpy as np
from ase.io import read

@pytest.fixture
def input_pystruc():
	py_struc = Structure.from_file("./data/Si.cif")
	return py_struc

@pytest.fixture
def input_pystruc2():
	py_struc2 = Structure.from_file("./data/LiF.cif")
	return py_struc2


@pytest.fixture
def ipt_LiF_forces(scope="session"):
	ase_struc = read("./data/LiF_p1.cif")
	atf = np.loadtxt("./data/LiF_p1_forces.txt")
	return ase_struc, atf


    
    
   