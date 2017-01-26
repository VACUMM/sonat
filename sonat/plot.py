"""Plot utilities"""
#
# Copyright IFREMER (2016-2017)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#

import os
import re
import matplotlib.pyplot as P
from matplotlib import rc_params_from_file, rcParams
from vcmq import (map, hov, stick, curve, section, dict_check_defaults, plot2d)

from .__init__ import SONATError, sonat_warn
from .misc import slice_gridded_var


#: Matplotlib default configuration file
SONAT_DEFAULT_MATPLOTLIBRC =  os.path.join(os.path.dirname(__file__), 'matplotlibrc')

#: Matplotlib user configuration file
SONAT_USER_MATPLOTLIBRC =  'matplotlibrc'

#: Default arguments to plots
DEFAULT_PLOT_KWARGS = dict(
    contour_linewidths=.5,
    colorbar_fraction=.1,
    quiver_units='dots',
    quiver_width=1.2,
    fill='contourf',
    quiver_norm=3,
    proj='merc',
    cmap='auto',
    autoresize='y',
    autoresize_minaspect=.5,
    colorbar_shrink=.8,
    fillcontinents_zorder=10,
#    figsize=(5, 3), # in matplotlibrc
    show=False,
    close=True,
    )

RE_GRIDDED_ORDER_MATCH = re.compile(r'\-?t?z?y?x?$').match

def plot_gridded_var(var, member=None, time=None, depth=None, lat=None, lon=None,
        plotfunc=None, **kwargs):
    """Generic 1D or 2D plot of a [T][Z]YX variable

    Parameters
    ----------
    var: MV2.array to tuple of that
        The gridded variable. If it has a vertical dimension, coordinates
        must be 1D.
    time: float, time, None, slice
        Interpolate var at this time
    depth: float, slice, None
        Make an horizontal section at this negative depth
    lat: float, slice, None
        Make a zonal section at this latitude
    lon: float, slice, None
        Make meridional section at this longitude
    """
    # Tuple
    vv = var[:2] if isinstance(var, tuple) else (var, )
    quiver = len(vv) == 2
    var = vv[0]

    # Check order
    order = var.getOrder()
    if not RE_GRIDDED_ORDER_MATCH(order):
        raise SONATError('Wrong order for variable: {}. '
            'It must match [-][t][z][y][x]'.format(order))

    # Slice
    vv = tuple([slice_gridded_var(var, member=member, time=time, depth=depth,
            lat=lat, lon=lon)(squeeze=1)  for var in vv])

    # Scalar
    for var in vv:
        if var.ndim==0:
            sonat_warn('Cannot plot a scalar')
            return

    # Make sure to have at most 2 dims
    vv = list(vv)
    for iv, var in enumerate(vv):
        if var.ndim>2:
            order = var.getOrder()
            for i in xrange(var.ndim-2):
                if order[0]=='z':
                    var = var[0]
                else:
                    var = var[-1]
            vv[iv] = var

    # Check coherence
    if not quiver:
        vv = var
    elif vv[0].shape != vv[1].shape:
        raise SONATError('Variables to plot have incompatible shapes')

    # Select plot function
    if plotfunc is None:
        order = var.getOrder()
        if var.ndim==1:
            plotfunc = stick if quiver else curve
        elif 't' in order:
            plotfunc = hov
        elif 'z' in order:
            plotfunc = section
        elif order == 'yx':
            plotfunc = map
        else:
            plotfunc = plot2d

    # Positional arguments
    if plotfunc is stick:
        args = vv
    else:
        args = [vv]

    # Default optional arguments
    dict_check_defaults(kwargs,  **DEFAULT_PLOT_KWARGS)

    # Plot and return
    return plotfunc(*args, **kwargs)


def load_mplrc(userfile=None):
    """Load a matplotlib or default user configuration file"""
    # Load default file first
    rcParams.update(rc_params_from_file(SONAT_DEFAULT_MATPLOTLIBRC, use_default_template=False))

    # Load user file
    userfile = str(userfile)
    if userfile=='False':
        return
    if userfile=='True':
        userfile = 'None'
    if userfile=='None':
        userfile = SONAT_USER_MATPLOTLIBRC
    if not os.path.exists(userfile):
        return
    rcParams.update(rc_params_from_file(userfile, use_default_template=False))

load_mplrc()
