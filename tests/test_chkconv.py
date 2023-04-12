# -*- coding: utf-8 -*-
"""Testing"""
# import numpy as np
# import pytest
from musConv.chkconv import ChkConvergence


def test_force_length_exception(ipt_lif_forces):
    """test force length"""
    try:
        ChkConvergence(ipt_lif_forces[0], ipt_lif_forces[1].pop(-1))
        assert False
    except ValueError:
        assert True


def test_mu_var_exception(ipt_lif_forces):
    """test mu variable"""
    try:
        ChkConvergence(ipt_lif_forces[0], ipt_lif_forces[1], "Fe")
        assert False
    except ValueError:
        assert True


def test_conv_cond(ipt_lif_forces):
    """test conv condition"""
    csc = ChkConvergence(ipt_lif_forces[0], ipt_lif_forces[1], "H")
    cond = csc.apply_first_crit()
    cond2 = csc.apply_2nd_crit()
    assert cond is False and all(cond2) is False
