import pytest
import numpy as np
from musConv.chkconv import chkSCconvergence


def test_force_length_exception(ipt_LiF_forces):
	try:
		chkSCconvergence(ipt_LiF_forces[0],ipt_LiF_forces[1].pop(-1))
		assert False
	except:
		assert True


def test_mu_var_exception(ipt_LiF_forces):
	try:
		chkSCconvergence(ipt_LiF_forces[0],ipt_LiF_forces[1], 'Fe')
		assert False
	except:
		assert True


def test_conv_cond(ipt_LiF_forces):
	csc   = chkSCconvergence(ipt_LiF_forces[0],ipt_LiF_forces[1], 'H')
	cond  = csc.apply_first_crit()
	cond2 = csc.apply_2nd_crit()
	assert cond == False and all(cond2) == False

