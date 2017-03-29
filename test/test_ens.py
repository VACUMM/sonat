"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4, create_dep, func_name

from util import (THISDIR, NCPAT_MANGA, assert_allclose, LOGGER, NCFILE_MANGA0,
    NCFILE_MANGA1)

from sonat.ens import (load_model_at_regular_dates, generate_pseudo_ensemble,
    Ensemble)

ENS_NCFILE = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')

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
    varnames = ['temp', 'depth']
    level = None
    temp, depth = load_model_at_regular_dates(ncpat, varnames=varnames, time=time,
        lat=lat, lon=lon, level=level, modeltype='mars', nt=nt, dtfile=dtfile, sort=True)
    assert temp.shape==(nt, 15, 10, 4)
    assert depth.shape==(nt, 15, 10, 4)
    ctimes = temp.getTime().asComponentTime()
    assert ctimes[1] == comptime('2014-1-3 0:0:0.0')
    assert ctimes[-2] == comptime('2014-1-24 0:0:0.0')

    # 3D with z interpolation
    depths = create_dep([-40., -30, -20, -10, 0.])
    temp = load_model_at_regular_dates(ncpat, varnames=varnames[0], time=time,
        lat=lat, lon=lon, modeltype='mars', nt=2, dtfile=dtfile, sort=True,
        depths=depths)
    assert temp.shape == (2, len(depths), 10, 4)

    # Surf
    varnames = 'temp'
    level = {'temp':'surf'}
    temp = load_model_at_regular_dates(ncpat, varnames=varnames, time=time,
        lat=lat, lon=lon, level=level, modeltype='mars', nt=nt, dtfile=dtfile, sort=True)
    assert temp.shape==(nt, 10, 4)
    assert temp.id == 'temp_surf'

def test_ens_generate_pseudo_ensemble():

    # Specs
    ncpat = NCPAT_MANGA
    varnames = ['temp', 'sal', 'u', 'v']
    time = ('2014-01-01 13', '2014-01-25 12')
    nrens = 14
#    nrens = 2
    enrich = 1.5
    dtfile = (12, 'day')
    ncfile = ENS_NCFILE
    level = {'temp':('3d', 'surf'), 'u':'surf', 'v':'surf'}
    depths = create_dep([-40., -30, -20, -10, 0.])

    # Direct
    enrich = 0 # <= 1
    (temp, temp_surf, sal, u_surf, v_surf) = generate_pseudo_ensemble(ncpat,
        nrens=nrens, enrich=enrich,
        time=time, varnames=varnames, dtfile=dtfile, logger=LOGGER, anomaly=False,
        level=level, depths=depths)
    assert temp.shape[0]==nrens
    assert sal.shape[0]==nrens
    assert temp.ndim==4
    assert temp.shape[1]==len(depths)
    assert temp_surf.ndim==3
    assert v_surf.ndim==3
    f = cdms2.open(NCFILE_MANGA0)
    temp0 = f('TEMP', time=slice(1, 2), level=slice(-1, None), squeeze=1)
    f.close()
    assert_allclose(temp0, temp_surf[0])
    tsum = temp.sum()

    # Enrichment
    enrich = 1.5
    ens = generate_pseudo_ensemble(ncpat, nrens=nrens, enrich=enrich,
        varnames=varnames, time=time, dtfile=dtfile, logger=LOGGER, getmodes=True,
        level=level, depths=depths)
    (temp, temp_surf, sal, u_surf, v_surf), modes = ens
    (temp_eof, temp_surf_eof, sal_eof, u_surf_eof, v_surf_eof) = modes['eofs']
    ev = modes['eigenvalues']
    temp_var, temp_surf_var, sal_var, u_surf_var, v_surf_var = modes['variance']
    assert tsum!=temp.sum()
    assert temp.shape[0]==nrens
    assert sal.shape[0]==nrens
    eof0 = N.concatenate( (temp_eof[0].compressed(), temp_surf_eof[0].compressed(),
        sal_eof[0].compressed(),
        u_surf_eof[0].compressed(), v_surf_eof[0].compressed()))
    assert_allclose((eof0**2).sum(), 1)
    eof1 = N.concatenate( (temp_eof[1].compressed(), temp_surf_eof[1].compressed(),
        sal_eof[1].compressed(),
        u_surf_eof[1].compressed(), v_surf_eof[1].compressed()))
    assert_allclose((eof0*eof1).sum(), 0, atol=1e-7)
    assert_allclose(ev.total_variance, eof0.size)
    expv = (ev**2).sum()/ev.total_variance
    assert expv > .8 and expv < 1
    expvm = temp.var(axis=0).mean()/temp_var.mean()
    assert expvm > .8 and expvm < 1

    # Save ensemble
    f = cdms2.open(ncfile, 'w')
    for var in (
            temp, temp_surf, sal, u_surf, v_surf,
            temp_eof, temp_surf_eof, sal_eof, u_surf_eof, v_surf_eof,
            temp_var, temp_surf_var, sal_var, u_surf_var, v_surf_var,
            ev
            ):
        f.write(var)
    f.close()

def test_ens_ensemble_init():

    # Load file from previous routine
    ncfile = ENS_NCFILE
    varnames = ['sal', 'temp', 'temp_surf', 'u_surf',  'v_surf']
    vars = []
    f = cdms2.open(ncfile)
    for vname in f.listvariables():
        if vname in varnames:
            vars.append(f(vname))
    f.close()
    varnames = [v.id for v in vars]

    # Init from variables
    ensv = Ensemble(vars, checkvars=True)

    # Init from file
    ensf = Ensemble.from_file(ncfile, checkvars=True, logger=LOGGER)

    # Checks
    assert [v.id for v in ensv.variables] == varnames
    assert [v.id for v in ensf.variables] == varnames
    assert_allclose(ensv.stacked_data, ensf.stacked_data)


CACHE = {}
def get_ens():

    # From cache
    if 'ens' in CACHE:
        return CACHE['ens']

     # Load from file
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')
    ens = Ensemble.from_file(ncfile, logger=LOGGER)
    CACHE['ens'] = ens
    return ens

def test_ens_ensemble_get_diags():

    # Get ens
    ens = get_ens()

    # Diags
    diags = ens.get_diags()

def test_ens_ensemble_plot_diags():

    # Get ens
    ens = get_ens()

    # Diags
    figs = ens.plot_diags(
        zonal_sections=[47.5], merid_sections=[-4],
        kurtosis=False, normaltest=False, skewtest=False,
        kurtosistest=False, skew=False, mean=False,
        )

def test_ens_ensemble_export_html_diags():

    # Get ens
    ens = get_ens()

    # Plot and export diags
    ens.export_html_diags(func_name()+'.html', surf=True, variance=True,
        mean=False, skewtest=False, kurtosistest=False, normaltest=False)



if __name__=='__main__':
    test_ens_load_model_at_regular_dates()
    test_ens_generate_pseudo_ensemble()
    test_ens_ensemble_init()
    test_ens_ensemble_get_diags()
    test_ens_ensemble_plot_diags()
    test_ens_ensemble_export_html_diags()
