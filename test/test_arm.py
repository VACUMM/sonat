"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4, P, func_name

from util import (THIS_DIR, NCPAT_MANGA, assert_allclose, LOGGER, NCFILE_MANGA0,
    NCFILE_MANGA1, NCFILE_OBS_HFRADARS, NCFILE_OBS_PROFILES, NCFILE_OBS_SATSST,
    get_bathy)

from sonat.ens import Ensemble
from sonat.obs import NcObsPlatform, ObsManager
from sonat.arm import (ARM, register_arm_score_function, get_arm_score_function,
    ARM_SCORE_FUNCTIONS, XYLocARMSA)

netcdf4()

_CACHE = {}

def get_arm():
    if 'arm' not in _CACHE:
        _CACHE['arm'] = test_arm_arm_init()
    return _CACHE['arm']

def test_arm_arm_init():

    # Load ensemble
    ncfile = os.path.join(THIS_DIR, 'test_ens_generate_pseudo_ensemble.nc')
    ens = Ensemble.from_file(ncfile, checkvars=True, logger=LOGGER)

    # Load observations
    obs0 = NcObsPlatform(NCFILE_OBS_HFRADARS)
    obs1 = NcObsPlatform(NCFILE_OBS_PROFILES)
    obs2 = NcObsPlatform(NCFILE_OBS_SATSST)
    obsmanager = ObsManager([obs0, obs1, obs2])

    # Bathymetry
    bathy = get_bathy()[::2, ::2]

    # Init ARM
    arm = ARM(ens, obsmanager, syncnorms=True, bathy=bathy)

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

    # Raw spectrum
    arm.plot_spect()

def test_arm_arm_plot_arm():

    # Load ARM
    arm = get_arm()

    # Raw array modes
    arm.plot_arm(modes=1)

def test_arm_arm_plot_rep():

    # Load ARM
    arm = get_arm()

    # Raw array modes
    arm.plot_rep(modes=1, surf=True, zonal_sections=47.5,
                 sync_vminmax=False, nmax=25,
                 obs_lat_interval_width=.3, obs_legend_loc='upper right')

def test_arm_arm_indirect_spectrum():

    # Load ARM
    arm = get_arm()

    # Indirect spectrum
    pcs = N.dot(arm.S.T, arm.raw_arm)
    spect = (pcs**2).sum(axis=0)

    # Compare with direct spectrum
    assert_allclose(spect, arm.raw_spect, atol=1e-7)


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


def test_arm_scores():

    # Fake data
    spect = N.array([4.5, 3.5, 2.5, 1.5, .5, 0.1])
    arm = None
    rep = None

    # Scores
    res = [get_arm_score_function(score_type)(spect, arm, rep)
           for score_type in ('nev', 'fnev', 'relvar')]

    return res

def test_arm_arm_export_html():

    # Load ARM
    arm = get_arm()

    # Raw array modes
    arm.export_html(func_name()+'.html', obs_legend_loc='upper right',
                    modes=1, varnames=['temp'],
                    arm_full2d=False)

def test_arm_xylocarmsa():

     # Load ARM
    arm = get_arm()

    # Init sensivity analyser
    armsa = XYLocARMSA(arm)

    # Sensitivity analysis
    resd = armsa.analyse(direct=True, score_type='relvar')
    resi = armsa.analyse(direct=False, score_type='relvar')
    pass

def test_arm_xylocarmsa_plot():

     # Load ARM
    arm = get_arm()

    # Init sensivity analyser
    armsa = XYLocARMSA(arm)

    # Plot
    armsa.plot(score_type='fnev', direct=True)
    armsa.plot(score_type='fnev', direct=False)
    armsa.plot(score_type='relvar', direct=True)
    armsa.plot(score_type='relvar', direct=False)

def test_arm_xylocarmsa_export_html():

     # Load ARM
    arm = get_arm()

    # Init sensivity analyser
    armsa = XYLocARMSA(arm)

     # Plot
    armsa.export_html(func_name()+'.html', score_type='fnev')


if __name__=='__main__':
    res = test_arm_arm_init()
    res = test_arm_arm_project_ens_on_obs()
    res = test_arm_arm_inputs()
    res = test_arm_arm_analyse()
    res = test_arm_arm_results()
    res = test_arm_arm_indirect_spectrum()
    res = test_arm_register_arm_score_function()
    res = test_arm_get_arm_score_function()
    res = test_arm_scores()
    res = test_arm_arm_plot_spect()
    res = test_arm_arm_plot_arm()
    res = test_arm_arm_plot_rep()
    res = test_arm_arm_export_html()
    res = test_arm_xylocarmsa()
    res = test_arm_xylocarmsa_plot()
    res = test_arm_xylocarmsa_export_html()



