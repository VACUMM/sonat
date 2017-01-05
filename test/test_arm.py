"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4

from util import (THISDIR, NCPAT_MANGA, assert_allclose, LOGGER, NCFILE_MANGA0,
    NCFILE_MANGA1, NCFILE_OBS_SURF)

from sonat.ens import Ensemble
from sonat.obs import NcObsPlatform, ObsManager
from sonat.arm import ARM

netcdf4()

def test_arm_arm_init():

    # Load ensemble
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')
    ens = Ensemble.from_file(ncfile, checkvars=True)

    # Load observations
    obs_surf0 = NcObsPlatform(NCFILE_OBS_SURF, lon=(-7, -5), varnames=['temp'],
        norms=0.3)
    obs_surf1 = NcObsPlatform(NCFILE_OBS_SURF, lon=(-5, 0), varnames=['temp', 'sal'],
        norms=[0.2, 0.1])
    obsmanager = ObsManager([obs_surf0, obs_surf1])

    # Init ARM
    arm = ARM(ens, obsmanager, syncnorms=True)

    # Checks
    assert_allclose(arm.obsmanager[0].norms + arm.obsmanager[1].norms[:1],
        arm.ens.norms[0])

    return arm

def test_arm_arm_project_ens_on_obs():

    # Load ARM
    arm = test_arm_arm_init()

    # Get matrices
    oens = arm.project_ens_on_obs()

    return oens

def test_arm_arm_get_matrices():

    # Load ARM
    arm = test_arm_arm_init()

    # Get matrices
    mats = arm.get_matrices()

    return mats

def test_arm_arm_analyse():

    # Load ARM
    arm = test_arm_arm_init()

    # Analyse
    ana = arm.analyse()

    return ana



if __name__=='__main__':
    res = test_arm_arm_init()
    res = test_arm_arm_project_ens_on_obs()
    res = test_arm_arm_get_matrices()
    res = test_arm_arm_analyse()



