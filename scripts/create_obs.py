#!/usr/bin/env python
"""Create fake observation netcdf files
"""

from vcmq import (cdms2, MV2, os, create_dep, create_lon, create_lat, N,
                  masked_polygon, map2)

cdms2.setAutoBounds(False)

# Profiles
ncfile = '../data/obs.profiles.nc'
lons = MV2.array([-5.8, -5.7, -4.6, -2.8], id='lon')
lats = MV2.array([48.1, 47.5, 47.4, 47.3], id='lat')
daxis = create_dep((-100., 1, 5.))
nz = len(daxis)
np = len(lons)
temp_error = N.resize([.2, .5, .3, .8], (nz, np))
sal_error = N.resize([.3, .1, .15, .4], (nz, np))
temp_error = MV2.array(temp_error, id='temp_error')
sal_error = MV2.array(sal_error, id='sal_error')
temp_error[:nz/2, 2] = MV2.masked
sal_error[:nz/2, 2] = MV2.masked
temp_error[:3*nz/4, 3] = MV2.masked
sal_error[:3*nz/4, 3] = MV2.masked
mob = MV2.array([0, 1, 1, 0], id='mobility', fill_value=-1)
paxis = lons.getAxis(0)
paxis.id = 'station'
axes = [daxis, paxis]
f = cdms2.open(ncfile , 'w')
for var in lons, lats, mob, temp_error, sal_error:
    var.setAxis(-1, paxis)
    if var.ndim==2:
        var.setAxis(0, daxis)
    f.write(var)
f.close()
print os.path.abspath(ncfile)


# HF radars
ncfile = '../data/obs.hfradars.nc'
lon = create_lon((-5.8, -4.7, .1))
lat = create_lat((47.8, 48.3, .1))
nx = len(lon)
ny = len(lat)
u_error = MV2.ones((ny, nx), id='u_error', axes=[lat, lon])
u_error[:] *= .2
u_error[:] += N.resize(N.linspace(0, .1, nx), (ny, nx))
u_error[:2, :2] = u_error[-2:, :2] = N.ma.masked
v_error = u_error.clone()
v_error.id = 'v_error'
mob = MV2.zeros((ny, nx), id='mobility', axes=[lat, lon], dtype='i')
f = cdms2.open(ncfile , 'w')
for var in mob, u_error, v_error:
    f.write(var)
f.depth = 'surf'
f.close()
print os.path.abspath(ncfile)


# Satellite SST
ncfile = '../data/obs.satsst.nc'
lon = create_lon((-7, -2.1, .5))
lat = create_lat((46., 49.1, .5))
nx = len(lon)
ny = len(lat)
sst_error = MV2.ones((ny, nx), id='temp_error', axes=[lat, lon])
sst_error[:] *= .5
sst_error = masked_polygon(sst_error, 'i')
f = cdms2.open(ncfile , 'w')
f.write(sst_error)
f.mobility = 0.
f.depth = 'surf'
f.close()
print os.path.abspath(ncfile)

