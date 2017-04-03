"""Test script for module :mod:`sonat.plot`"""

import os
import cdms2
#from matplotlib import rcParams
from vcmq import adatetime, comptime, func_name, P

from util import (NCFILE_MANGA0, THIS_DIR)

from sonat.plot import (plot_gridded_var, create_map)

def test_plot_gridded_var():
    # Get var
    ncfile = os.path.join(THIS_DIR, 'test_ens_generate_pseudo_ensemble.nc')
    f = cdms2.open(ncfile)
    temp = f('temp')
    temp_surf = f('temp_surf')
    sal = f('sal')
    u = f('u_surf')
    v = f('v_surf')
    f.close()

    # Figure file
    figpat = os.path.join(THIS_DIR, func_name()) + '.{}.png'

    # Plots
    kw = dict(show=False)
   # - 1D
    kw1d = dict(kw, tight_layout=True, figsize=(6, 3))
    plot_gridded_var(temp_surf, lat=47.5, lon=-4.,
        savefig=figpat.format('m.scalar'), **kw1d)
    plot_gridded_var(temp_surf, member=slice(0, 1), lon=-4.,
        savefig=figpat.format('y.scalar'), **kw1d)
    plot_gridded_var(temp_surf, lat=47.5, member=1,
        savefig=figpat.format('x.scalar'), **kw1d)
    plot_gridded_var((u, v), lat=47.5, member=1, depth=0,
        savefig=figpat.format('x.tuple'), **kw1d)
    # - 2D
    kw2d = dict(kw, figsize=(6, 4))
    plot_gridded_var(temp_surf, member=slice(0, 1), levels_mode='symetric',
        savefig=figpat.format('xy.scalar'), cmap='balance', **kw2d)
    plot_gridded_var(sal, member=slice(0, 1), depth=-20, levels_mode='symetric',
        savefig=figpat.format('xyh.scalar'), cmap='delta', **kw2d)
    plot_gridded_var((u, v), member=slice(0, 1),
        savefig=figpat.format('xy.tuple'), colorbar=False, fill=False,
        contour=False, title='Currents', cmap='speed', **kw2d)
    plot_gridded_var(temp, lat=47.5, levels_mode='symetric',
        savefig=figpat.format('mx.scalar'), cmap='balance', **kw2d)


def test_plot_create_map():

    # Read bathy
    f = cdms2.open(NCFILE_MANGA0)
    bathy = -f('H0')
    f.close()
    alon = bathy.getLongitude()
    alat = bathy.getLatitude()
    lon = (alon[:].min(), alon[:].max())
    lat = (alat[:].min(), alat[:].max())

    # 2D maps
    # - simple
    m = create_map(lon, lat)
    P.savefig(func_name()+'.2d.simple.png', title="2D / simple")
    P.close()
    # - bathy
    m = create_map(lon, lat, bathy=bathy, title="2D / bathy")
    P.savefig(func_name()+'.2d.bathy.png')
    P.close()

    # 3D maps
    # - simple
    m = create_map(lon, lat, level=-200, title="3D / simple")
    P.savefig(func_name()+'.3d.simple.png')
    P.close()
    # - bathy
    m = create_map(lon, lat, axes="3d", bathy=bathy, level=-200,
                   title="3D / bathy")
    P.savefig(func_name()+'.3d.bathy.png')
    P.close()


if __name__=='__main__':
    test_plot_gridded_var()
    test_plot_create_map()
