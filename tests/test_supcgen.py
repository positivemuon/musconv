# -*- coding: utf-8 -*-
"""Testing"""
import numpy as np
import pytest

from musconv.supcgen import ScGenerators


def test_min_length_except(input_pystruc):
    """Test min length"""
    scg = ScGenerators(input_pystruc)
    try:
        scg.initialize(np.min(input_pystruc.lattice.abc) - 1)
        assert False
    except ValueError:
        assert True


@pytest.mark.parametrize("input_struc", [("input_pystruc"), ("input_pystruc2")])
def test_gen_nearcubic_supc(input_struc, request):
    """supecell generator"""
    input_struc = request.getfixturevalue(input_struc)
    scg = ScGenerators(input_struc)
    py_sup, sc_matrix = scg.gen_nearcubic_supc(
        input_struc,
        input_struc.num_sites + 1,
        np.Inf,
        np.min(input_struc.lattice.abc) + 1,
    )
    assert np.linalg.det(sc_matrix) != 0
    assert np.min(py_sup.lattice.abc) > np.min(input_struc.lattice.abc)
    assert np.min(py_sup.lattice.angles) >= np.min(input_struc.lattice.angles)


@pytest.mark.parametrize("input_struc", [("input_pystruc"), ("input_pystruc2")])
def test_gen_larger_supc(input_struc, request):
    """Test larger supercell generator"""
    input_struc = request.getfixturevalue(input_struc)
    scg = ScGenerators(input_struc)
    st_m, sc_m, mu_c = scg.initialize()
    st2_m, sc2_m = scg.re_initialize(st_m, mu_c)
    assert st2_m.num_sites > st_m.num_sites
    assert np.linalg.det(sc2_m) > np.linalg.det(sc_m)
