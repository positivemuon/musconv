# -*- coding: utf-8 -*-
""" Run example"""
from pymatgen.core import Structure

from musconv.supcgen import ScGenerators

if __name__ == "__main__":
    # load structure with pymatgen
    py_struc = Structure.from_file("LiF.cif")

    sg = ScGenerators(py_struc)

    # initialize the caluclations
    # py_scst_mu2,sc_matrix,mu_frac_coord=sg.initialize(min_length)
    py_scst_mu2, sc_matrix, mu_frac_coord = sg.initialize()
    py_scst_mu2.to(filename="positions.cif".format())
    # print(sc_matrix)

    # while and if loop then depending on workchain usage
    # py_scst_mu2,sc_matrix=sg.re_initialize(py_scst_mu2,mu_frac_coord)
