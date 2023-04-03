import pytest
import numpy as np
from musConv.supcgen import SCgenerators


def test_min_length_except(input_pystruc):
	sc = SCgenerators(input_pystruc)
	try:
		sc.initialize(np.min(input_pystruc.lattice.abc)-1)
		assert False
	except:
		assert True


@pytest.mark.parametrize("input_struc", [("input_pystruc"), ("input_pystruc2")])
def test_gen_nearcubic_SC(input_struc,request):
	input_struc = request.getfixturevalue(input_struc)
	sc = SCgenerators(input_struc)
	py_sup, sc_matrix = sc.gen_nearcubic_SC(
		input_struc, 
		input_struc.num_sites+1, 
		np.Inf, 
		np.min(input_struc.lattice.abc)+1
		)
	assert np.min(py_sup.lattice.abc) > np.min(input_struc.lattice.abc)
	assert np.min(py_sup.lattice.angles) >= np.min(input_struc.lattice.angles)


@pytest.mark.parametrize("input_struc", [("input_pystruc"), ("input_pystruc2")])
def test_gen_larger_supc(input_struc,request):
	input_struc = request.getfixturevalue(input_struc)
	sc = SCgenerators(input_struc)
	st_m, sc_m, mu_c   = sc.initialize()
	st2_m, sc2_m   = sc.re_initialize(st_m, mu_c)
	assert st2_m.num_sites > st_m.num_sites
	assert np.linalg.det(sc2_m) > np.linalg.det(sc_m)
		
