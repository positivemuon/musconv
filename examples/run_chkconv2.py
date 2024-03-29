# -*- coding: utf-8 -*-
""" Run example"""
import random
import sys

import numpy as np
from ase.atom import Atom
from ase.calculators.lj import LennardJones
from ase.data import atomic_numbers
from ase.io import read

from musconv.chkconv import ChkConvergence

try:
    atoms = read(sys.argv[1])
except:
    sys.exit()

valid_atom = set(atomic_numbers.keys()) - set(atoms.get_chemical_symbols())

IMP = random.sample(valid_atom, 1)[0]

atoms.append(Atom(IMP, np.random.random(3)))
atoms.set_calculator(LennardJones())

csc = ChkConvergence(atoms, atoms.get_forces(), mu_num_spec=IMP)
COND = csc.apply_first_crit()
cond2 = csc.apply_2nd_crit()
print(f"Convergence of 1st criteria is {COND}, while 2nd criteria is {cond2}")
