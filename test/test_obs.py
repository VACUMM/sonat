"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4, DS

from util import (assert_allclose, LOGGER, NCFILE_OBS_SURF, NCFILE_MANGA0,
    NCFILE_OBS_SURF_UV, NCFILE_OBS_3D_TS)

from sonat.obs import (NcObsPlatform, ObsManager, load_obs,
    register_obs_platform, get_obs_platform, OBS_PLATFORM_TYPES,
    load_obs_platform)


def test_obs_ncobsplatform_surf():

    # Load and stack obs
    obs = NcObsPlatform(NCFILE_OBS_SURF, logger=LOGGER, lat=(45, 47.8),
        norms=[0.2, 0.1])
    stacked = obs.stacked_data.copy()
    assert obs.lons.shape == (3, )
    assert obs.ns == 6
    assert obs.ndim == 1
    assert_allclose(obs.means, [0, 0], atol=1e-7)
    assert obs.depths == 'surf'
    assert_allclose(stacked, [ 2.5,  1.5,  4. ,  1. ,  1.5,  4. ])

    # Named norms
    notnormed = obs.set_named_norms(temp=0.1)
    assert obs.norms == [0.1, 0.1]
    assert_allclose(obs.stacked_data[:3], 2*stacked[:3])
    assert_allclose(obs.stacked_data[3:], stacked[3:])
    obs.set_named_norms(temp=0.2)

    # Interpolate model
    f = DS(NCFILE_MANGA0, 'mars', level=obs.depths)
    temp = f('temp')
    sal = f('sal')
    f.close()
    otemp = obs.project_model(temp)
    osal = obs.project_model(sal)
    otem_true = [12.97311556, 12.91558515, 10.58179214]
    assert_allclose(otemp[0], otem_true)

    # Stack model
    otemp[:] -= 11.5
    osal[:] -= 35.5
    stacked_data = obs.restack([otemp, osal])
    assert stacked_data.shape == (6, 15)
    assert_allclose(stacked_data[:3,0]*obs.norms[0] + 11.5, otem_true)

    return obs


def test_obs_ncobsplatform_surf_gridded():

    # Load and stack obs
    obs = NcObsPlatform(NCFILE_OBS_SURF_UV)
    stacked = obs.stacked_data.copy()
    assert stacked.ndim==1


def test_obs_obsmanager_init():

    # Load and stack surface obs
    obs_surf0 = NcObsPlatform(NCFILE_OBS_SURF, lon=(-7, -5), varnames=['temp'],
        norms=0.2)
    obs_surf1 = NcObsPlatform(NCFILE_OBS_SURF, lon=(-5, 0), varnames=['temp', 'sal'],
        norms=[0.2, 0.1])

    # Setup manager
    manager = ObsManager([obs_surf0, obs_surf1])

    # Checks
    stacked = manager.stacked_data
    assert_allclose(stacked, [ 1. ,  2.5,  1.5,  4. ,  1.5,  4. ])
    assert_allclose(manager.lons, [-5.8,  -5.7,  -4.6,  -2.8])
    assert_allclose(manager.lats, [48.1, 47.5, 47.4, 47.3])
    assert manager.varnames == ['temp', 'sal']
    model_specs = manager.get_model_specs()
    assert sorted(model_specs.keys()) == ['depths', 'lat', 'lon', 'varnames']
    assert model_specs['varnames'] == ['temp', 'sal']
    assert model_specs['depths']['temp'] == ('surf', )
    assert model_specs['lat'][:2] == (47.3, 48.1)
    assert model_specs['lon'][:2] == (-5.8, -2.8)

    # Renorm by name
    manager.set_named_norms(temp=0.1)
    assert manager.stacked_data[0] == 2 * stacked[0]

    return manager

def test_obs_load_obs():

    # Setup manager
    manager = load_obs([NCFILE_OBS_3D_TS, NCFILE_OBS_SURF_UV])

    return manager

def test_obs_obsmanager_project_model():

    # Load manager
    manager = test_obs_obsmanager_init()

    # Interpolate model
    f = DS(NCFILE_MANGA0, 'mars', level=manager.obsplats[0].depths)
    temp = f('temp')
    f.close()
    otemp = manager.project_model(temp)
    assert_allclose(otemp[1][0], [12.91558515, 10.58179214])

    return manager


def test_obs_register_obs_platform():

    class MyNcObsPlatform(NcObsPlatform):
        platform_type = 'myplatform'

    register_obs_platform(MyNcObsPlatform, warn=False)

    assert MyNcObsPlatform.platform_type in OBS_PLATFORM_TYPES
    return MyNcObsPlatform

def test_obs_get_obs_platform():

    clso = test_obs_register_obs_platform()

    clsg = get_obs_platform('myplatform')

    assert clso is clsg


def test_obs_load_obs_platform():

    clso = test_obs_register_obs_platform()

    obs = load_obs_platform('myplatform', NCFILE_OBS_SURF)
    assert obs.varnames == ['temp', 'sal']
    return obs


if __name__=='__main__':
    res = test_obs_ncobsplatform_surf()
    res = test_obs_ncobsplatform_surf_gridded()
    res = test_obs_obsmanager_init()
    res = test_obs_load_obs()
    res = test_obs_obsmanager_project_model()
    res = test_obs_register_obs_platform()
    res = test_obs_get_obs_platform()
    res = test_obs_load_obs_platform()
