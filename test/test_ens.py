"""Test script for module :mod:`pyarm.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4

from util import (THISDIR, NCPAT_MANGA, assert_allclose, LOGGER, NCFILE_MANGA0,
    NCFILE_MANGA1)

from pyarm.ens import (load_model_at_regular_dates, generate_pseudo_ensemble,
    Ensemble)

netcdf4()

def test_ens_load_model_at_regular_dates():

    # Specs
    ncpat = NCPAT_MANGA
    time = ('2014-01-01 12', '2014-01-25 12')
    nt = 15
    lat = (47.2, 47.9)
    lon = (-5.5, -5)
    dtfile = (12, 'day')

    # 3D
    ncvars = ['temp', 'depth']
    level = None
    temp, depth = load_model_at_regular_dates(ncpat, ncvars=ncvars, time=time,
        lat=lat, lon=lon, level=level, modeltype='mars', nt=nt, dtfile=dtfile, sort=True)
    assert temp.shape==(nt, 15, 10, 4)
    assert depth.shape==(nt, 15, 10, 4)
    ctimes = temp.getTime().asComponentTime()
    assert ctimes[1] == comptime('2014-1-3 0:0:0.0')
    assert ctimes[-2] == comptime('2014-1-24 0:0:0.0')

    # Surf
    ncvars = 'temp'
    level = {'temp':'surf'}
    temp = load_model_at_regular_dates(ncpat, ncvars=ncvars, time=time,
        lat=lat, lon=lon, level=level, modeltype='mars', nt=nt, dtfile=dtfile, sort=True)
    assert temp.shape==(nt, 10, 4)

def test_ens_generate_pseudo_ensemble():

    # Specs
    ncpat = NCPAT_MANGA
    ncvars = ['temp', 'sal']
    time = ('2014-01-01 13', '2014-01-25 12')
    nrens = 14
    enrich = 1.5
    dtfile = (12, 'day')
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')

    # Direct
    enrich = 0 # <= 1
    (temp, sal) = generate_pseudo_ensemble(ncpat, nrens=nrens, enrich=enrich,
        time=time, ncvars=ncvars, dtfile=dtfile, logger=LOGGER, anomaly=False)
    assert temp.shape[0]==nrens
    assert sal.shape[0]==nrens
    f = cdms2.open(NCFILE_MANGA0)
    temp0 = f('TEMP', time=slice(1, 2), squeeze=1)
    f.close()
    assert_allclose(temp0, temp[0])

    # Enrichment
    enrich = 1.5
    ens = generate_pseudo_ensemble(ncpat, nrens=nrens, enrich=enrich,
        ncvars=ncvars, time=time, dtfile=dtfile, logger=LOGGER, getmodes=True)
    (temp, sal), modes = ens
    (temp_eof, sal_eof) = modes['eofs']
    ev = modes['eigenvalues']
    temp_var, sal_var = modes['variance']
    assert temp.shape[0]==nrens
    assert sal.shape[0]==nrens
    eof0 = N.concatenate( (temp_eof[0].compressed(), sal_eof[0].compressed()))
    assert_allclose((eof0**2).sum(), 1)
    eof1 = N.concatenate( (temp_eof[1].compressed(), sal_eof[1].compressed()))
    assert_allclose((eof0*eof1).sum(), 0, atol=1e-7)
    assert_allclose(ev.total_variance, eof0.size)
    expv = (ev**2).sum()/ev.total_variance
    assert expv > .8 and expv < 1
    expvm = temp.var(axis=0).mean()/temp_var.mean()
    assert expvm > .8 and expvm < 1

    # Save ensemble
    f = cdms2.open(ncfile, 'w')
    for var in temp, sal, temp_eof, sal_eof, ev, temp_var, sal_var:
        f.write(var)
    f.close()

def test_ens_ensemble_init():

    # Load file from previous routine
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')
    f = cdms2.open(ncfile)
    temp = f('temp')
    sal = f('sal')
    f.close()

    # Init from variables
    ensv = Ensemble([temp, sal])

    # Init from file
    ensf = Ensemble.from_file(ncfile)

    # Checks
    assert_allclose(ensv.stacked_data, ensf.stacked_data)

def test_ens_ensemble_diags():

    # Load from file
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')
    ens = Ensemble.from_file(ncfile)

    # Diags
    diags = ens.get_diags()
    pass


if __name__=='__main__':
#    test_ens_load_model_at_regular_dates()
    test_ens_generate_pseudo_ensemble()
#    test_ens_ensemble_init()
    test_ens_ensemble_diags()
