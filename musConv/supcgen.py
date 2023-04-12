# -*- coding: utf-8 -*-
"""Generates a nearly cubic supercell (SC) for convergence checks.
 DEPENDENCIES:
 impurity generator now in a  pymatgen extension
 (i) pymatgen-analysis-defects (pymatgen>=2022.10.22)
 (ii) numpy
"""
import numpy as np
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation,
)


class ScGenerators:
    """
    Generates a nearly cubic supercell (SC) for convergence checks.
    Inserts a muon in the supercell at a Voronoi interstitial site.
    One of it methods initializes the supercell generation and the other
    re-initializes generation of a larger supercell-size than the former.

    Param:
        py_struc: A pymatgen "unitcell" structure instance
                  This is used to create the supercell.
    """

    @staticmethod
    def gen_nearcubic_supc(py_struc, min_atoms, max_atoms, min_length):
        """
        Function that generates the nearly cubic supercell (SC).

        Params:
            py_struc         : The pymatgen structure instance
            min_atoms        : Integer-->Min number of atoms allowed in SC
            max_atoms        : Integer-->Max number of atoms allowed in SC
            min_length       : Integer-->Min length of the smallest SC lattice vector

        Returns:
            A nearly cubic SC structure and an array of the SC grid size
        """

        cst = CubicSupercellTransformation(
            min_atoms=min_atoms,
            max_atoms=max_atoms,
            min_length=min_length,
            force_diagonal=False,
        )

        py_scst = cst.apply_transformation(py_struc)
        sc_mat = cst.transformation_matrix

        return py_scst, sc_mat

    @staticmethod
    def append_muon_to_supc(py_scst, sc_mat, mu_frac_coord):
        """
        Add the muon as a hydrogen atom to the supercell (SC).

        Params:
            py_scst    : The pymatgen supercell structure
            sc_mat          : array-->the SC grid size
            mu_frac_coord     : array-->Interstitial site scaled in units

        Returns:
            A Pymatgen supercell structure that has the muon as a H atom at a Voronoi site


        """

        mu_frac_coord_sc = (np.dot(mu_frac_coord, np.linalg.inv(sc_mat))) % 1
        py_scst_withmu = py_scst.copy()

        # what if a H specie is in the structure object?
        try:
            py_scst_withmu.append(
                species="H",
                coords=mu_frac_coord_sc,
                coords_are_cartesian=False,
                validate_proximity=True,
            )
        except ValueError:
            raise SystemExit(
                "ValueError:The muon is too close to an existing site!, "
                "change muon site. Exiting...."
            ) from None

        return py_scst_withmu

    def __init__(self, py_struc):
        self.py_struc = py_struc
        self.max_atoms = np.Inf
        # self.py_scst     = None
        # self.mu_frac_coord  = None

    def initialize(self, min_length: float = None):
        """
        This func initializes the first supercell (SC) generation
        with the following conditions;

        min_atoms  : number of atoms in the given struc + 1
        max_atoms  : number of atoms in the given struc*(2**3)
                    This limits the SC generation to 8 times of the given cell.
        min_length : Min length of the smallest SC lattice vector

        Returns:
            A Pymatgen supercell structure that has the muon as a H atom at a Voronoi site
        """

        min_atoms = self.py_struc.num_sites + 1
        min_length = min_length or np.min(self.py_struc.lattice.abc) + 1

        if min_length < np.min(self.py_struc.lattice.abc):
            raise ValueError(
                " Provided supercell min_length is less than the length of the"
                " smallest input cell lattice vector"
            )

        py_scst, sc_mat = self.gen_nearcubic_supc(
            self.py_struc, min_atoms, self.max_atoms, min_length
        )

        # get a Voronoi interstitial site for the muon impurity
        vig = VoronoiInterstitialGenerator()
        mu_frac_coord = list(vig._get_candidate_sites(self.py_struc))[0][0]

        # Added 0.001 to move the impurity site from symmetric position
        mu_frac_coord = [x + 0.001 for x in mu_frac_coord]

        py_scst_with_mu = self.append_muon_to_supc(py_scst, sc_mat, mu_frac_coord)

        return py_scst_with_mu, sc_mat, mu_frac_coord

    def re_initialize(self, py_scst_with_mu, mu_frac_coord):
        """
        This function re-initializes the generation of a larger supercell-size in a loop
        when a condition is not met after the first initialization above.

        Param:
            iter_num : Integer--> iteration number in the loop

        Returns:
            A Pymatgen supercell structure that has the muon as a H atom at a Voronoi site
        """

        min_atoms = py_scst_with_mu.num_sites + 1
        min_length = np.min(py_scst_with_mu.lattice.abc) + 1

        py_scst, sc_mat = self.gen_nearcubic_supc(
            self.py_struc, min_atoms, self.max_atoms, min_length
        )

        py_scst_with_mu = self.append_muon_to_supc(py_scst, sc_mat, mu_frac_coord)

        return py_scst_with_mu, sc_mat


if __name__ == "__main__":
    import argparse

    from pymatgen.core import Structure

    parser = argparse.ArgumentParser(description="Generate nearly cubic supercell")
    parser.add_argument(
        "--min_length",
        metavar="N",
        type=float,
        default=None,
        help="Min length of the smallest SC lattice vector",
    )
    parser.add_argument("input_structure")

    args = parser.parse_args()

    min_len = args.min_length

    # load structure with pymatgen
    py_st = Structure.from_file(args.input_structure)

    sg = ScGenerators(py_st)

    # initialize the caluclations
    # py_scst_mu2,sc_mat,mu_frac_cd=sg.initialize(min_length)
    py_scst_mu2, sc_mt, mu_frac_cd = sg.initialize()
    py_scst_mu2.to(filename="positions.cif".format())
    # print(sc_mt)

    # while and if loop then depending on workchain usage
    # py_scst_mu2,sc_mt=sg.re_initialize(py_scst_mu2,mu_frac_cd)
