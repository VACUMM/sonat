"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4, DS

from util import (assert_allclose, LOGGER, NCFILE_OBS_SURF, NCFILE_MANGA0)

from sonat.obs import (NcObsPlatform, ObsManager)


def test_obs_ncobsplatform_surf():

    # Load and stack obs
    nop = NcObsPlatform(NCFILE_OBS_SURF, logger=LOGGER, lat=(45, 47.8),
        norms=[0.2, 0.1])
    assert nop.lons.shape == (3, )
    assert nop.ns == 6
    assert nop.ndim == 1
    assert_allclose(nop.means, [0, 0], atol=1e-7)
    assert nop.depths == 'surf'
    assert_allclose(nop.stacked_data, [ 2.5,  1.5,  4. ,  1. ,  1.5,  4. ])

    # Interpolate model
    f = DS(NCFILE_MANGA0, 'mars', level=nop.depths)
    temp = f('temp')
    sal = f('sal')
    f.close()
    otemp = nop.interp_model(temp)
    osal = nop.interp_model(sal)
    assert_allclose(otemp[0], [12.97311556, 12.91558515, 10.58179214])

    # Stack model
    otemp[:] -= 11.5
    osal[:] -= 35.5
    stacked_data = nop.restack([otemp, osal])

    return nop

def test_obs_obsmanager():

    # Load and stack surface obs
    obs_surf0 = NcObsPlatform(NCFILE_OBS_SURF, lon=(-7, -5), varnames=['temp'],
        norms=0.2)
    obs_surf1 = NcObsPlatform(NCFILE_OBS_SURF, lon=(-5, 0), varnames=['temp', 'sal'],
        norms=[0.2, 0.1])

    # Setup manager
    manager = ObsManager([obs_surf0, obs_surf1])
    assert_allclose(manager.stacked_data, [ 1. ,  2.5,  1.5,  4. ,  1.5,  4. ])
    assert_allclose(manager.lons, [-5.8,  -5.7,  -4.6,  -2.8])
    assert_allclose(manager.lats, [48.1, 47.5, 47.4, 47.3])
    assert manager.varnames == ['temp', 'sal']
    model_specs = manager.get_model_specs()

    # Interpolate model
    f = DS(NCFILE_MANGA0, 'mars', level=nop.depths)
    temp = f('temp')
    f.close()
    otemp = nop.interp_model(temp)

    return manager

if __name__=='__main__':
    res = test_obs_ncobsplatform_surf()
#    res = test_obs_obsmanager()
