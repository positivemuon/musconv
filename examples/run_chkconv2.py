# coding: utf-8
import sys
import random
from ase.io import read
from ase.data import atomic_numbers
from ase.atom import Atom
from ase.calculators.lj import LennardJones
import numpy as np
from musConv.chkconv import chkSCconvergence

try:
    atoms = read(sys.argv[1])
except:
    sys.exit()

valid_atom = set(atomic_numbers.keys()) - set(atoms.get_chemical_symbols())

imp = random.sample(valid_atom, 1)[0]

atoms.append(Atom(imp, np.random.random(3)))
atoms.set_calculator(LennardJones())

csc   = chkSCconvergence(atoms,atoms.get_forces(),mu_num_spec=imp)
cond  = csc.apply_first_crit()
cond2 = csc.apply_2nd_crit()
print(f"Convergence of 1st criteria is {cond}, while 2nd criteria is {cond2}")