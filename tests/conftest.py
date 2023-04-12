# -*- coding: utf-8 -*-
"""Testing"""
import numpy as np
import pytest
from ase.io import read
from pymatgen.core import Structure


@pytest.fixture
def input_pystruc():
    """input"""
    py_struc = Structure.from_file("./data/Si.cif")
    return py_struc


@pytest.fixture
def input_pystruc2():
    """input"""
    py_struc2 = Structure.from_file("./data/LiF.cif")
    return py_struc2


@pytest.fixture
def ipt_lif_forces(scope="session"):
    """input forces"""
    ase_struc = read("./data/LiF_p1.cif")
    atf = np.loadtxt("./data/LiF_p1_forces.txt")
    return ase_struc, atf
