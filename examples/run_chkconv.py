# -*- coding: utf-8 -*-
""" Run example"""
import numpy as np
from ase.io import read

from musconv.chkconv import ChkConvergence

if __name__ == "__main__":
    # load structure with ase and then forces
    ase_struc = read("LiF_p1.cif")
    atf = np.loadtxt("LiF_p1_forces.txt")

    # call the func
    csc = ChkConvergence(ase_struc, atf)
    COND = csc.apply_first_crit()
    cond2 = csc.apply_2nd_crit()
    print(f"Convergence of 1st criteria is {COND}, while 2nd criteria is {cond2}")
