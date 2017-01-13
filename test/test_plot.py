"""Test script for module :mod:`sonat.plot`"""

import os
import cdms2
#from matplotlib import rcParams
from vcmq import adatetime, comptime, func_name

from util import (NCFILE_MANGA0, THISDIR)

from sonat.plot import (plot_gridded_var, )

def test_plot_gridded_var():
    # Get var
    ncfile = os.path.join(THISDIR, 'test_ens_generate_pseudo_ensemble.nc')
    f = cdms2.open(ncfile)
    temp = f('temp')
    temp_surf = f('temp_surf')
    sal = f('sal')
    u = f('u_surf')
    v = f('v_surf')
    f.close()

    # Figure file
    figpat = os.path.join(THISDIR, func_name()) + '.{}.png'

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




if __name__=='__main__':
    test_plot_gridded_var()
