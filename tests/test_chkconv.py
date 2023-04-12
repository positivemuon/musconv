# -*- coding: utf-8 -*-
import numpy as np
import pytest

from musConv.chkconv import chkconvergence


def test_force_length_exception(ipt_lif_forces):
    try:
        chkconvergence(ipt_lif_forces[0], ipt_lif_forces[1].pop(-1))
        assert False
    except:
        assert True


def test_mu_var_exception(ipt_lif_forces):
    try:
        chkconvergence(ipt_lif_forces[0], ipt_lif_forces[1], "Fe")
        assert False
    except:
        assert True


def test_conv_cond(ipt_lif_forces):
    csc = chkconvergence(ipt_lif_forces[0], ipt_lif_forces[1], "H")
    cond = csc.apply_first_crit()
    cond2 = csc.apply_2nd_crit()
    assert cond == False and all(cond2) == False
