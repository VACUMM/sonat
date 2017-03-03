"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4, P

from util import (THISDIR, NCPAT_MANGA, assert_allclose, LOGGER, NCFILE_MANGA0,
    NCFILE_MANGA1, NCFILE_OBS_HFRADARS, NCFILE_OBS_PROFILES, NCFILE_OBS_SATSST)

from sonat.ens import Ensemble
from sonat.obs import NcObsPlatform, ObsManager
from sonat.arm import (ARM, register_arm_score_function, get_arm_score_function,
    ARM_SCORE_FUNCTIONS)

netcdf4()

_CACHE = {}

def get_arm():
    if 'arm' not in _CACHE:
        _CACHE['arm'] = test_arm_arm_init()
    return _CACHE['arm']

def test_arm_arm_init():

    # Load ensemble
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')
    ens = Ensemble.from_file(ncfile, checkvars=True, logger=LOGGER)

    # Load observations
    obs0 = NcObsPlatform(NCFILE_OBS_HFRADARS)
    obs1 = NcObsPlatform(NCFILE_OBS_PROFILES)
    obs2 = NcObsPlatform(NCFILE_OBS_SATSST)
    obsmanager = ObsManager([obs0, obs1, obs2])

    # Init ARM
    arm = ARM(ens, obsmanager, syncnorms=True)

    # Checks
    assert_allclose(arm.obsmanager[2].norms + arm.obsmanager[1].norms[:1],
        arm.ens.norms[arm.ens.varnames.index('temp')])

    return arm

def test_arm_arm_project_ens_on_obs():

    # Load ARM
    arm = get_arm()

    # Project
    oens = arm.project_ens_on_obs()

    return oens

def test_arm_arm_inputs():

    # Load ARM
    arm = get_arm()

    # Get matrices
    Yf = arm.Yf
    Af = arm.Af
    R = arm.R

    return

def test_arm_arm_analyse():

    # Load ARM
    arm = get_arm()

    # Analyse
    arm.analyse()

def test_arm_arm_results():

    # Load ARM
    arm = get_arm()

    # Raw results
    assert arm.raw_spect.shape == (arm.ndof, )
    assert arm.raw_arm.shape == (arm.nobs, arm.ndof)
    assert arm.raw_rep.shape == (arm.nstate, arm.ndof)


    return

def test_arm_arm_plot_spect():

    # Load ARM
    arm = get_arm()

    # Raw results
    arm.plot_spect()


def test_arm_register_arm_score_function():

    def arm_score_myfunc(ev, arm, rep):
        return 1

    register_arm_score_function(arm_score_myfunc)

    assert "myfunc" in ARM_SCORE_FUNCTIONS
    return arm_score_myfunc

def test_arm_get_arm_score_function():

    myfunc = test_arm_register_arm_score_function()

    func = get_arm_score_function('myfunc')

    assert func is myfunc


if __name__=='__main__':
#    res = test_arm_arm_init()
#    res = test_arm_arm_project_ens_on_obs()
#    res = test_arm_arm_inputs()
#    res = test_arm_arm_analyse()
#    res = test_arm_arm_results()
#    res = test_arm_register_arm_score_function()
#    res = test_arm_get_arm_score_function()
    test_arm_arm_plot_spect()



