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
from six import string_types
import numpy as N
import matplotlib.pyplot as P
from mpl_toolkits.mplot3d import art3d, Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import rc_params_from_file, rcParams
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize
import cdms2
from vcmq import (map, hov, stick, curve, section, dict_check_defaults, plot2d,
                  kwfilter, curv2rect, deplab, isaxis, meshbounds,
                  grid2xy, scalebox, Plot, land_color)

from .__init__ import SONATError, sonat_warn
from .misc import (slice_gridded_var, vminmax, slice_scattered_locs,
                   rescale_itv, get_long_name)


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
    contour_zorder=1,
    quiver_norm=3,
    proj='merc',
    cmap='auto',
    autoresize='y',
    autoresize_minaspect=.5,
    colorbar_shrink=.8,
#    fillcontinents_zorder=10,
#    figsize=(5, 3), # in matplotlibrc
    show=False,
    close=True,
    )

RE_GRIDDED_ORDER_MATCH = re.compile(r'\-?t?z?y?x?$').match

def plot_gridded_var(var, member=None, time=None, depth=None, lat=None, lon=None,
        plot_func=None, register_sm=True, **kwargs):
    """Generic 1D or 2D plot of a [T][Z]YX variable gridded variables

    Variables are sliced before plotting, using the
    :func:`~sonat.misc.slice_gridded_var` function.

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
    plot_func: vacumm plot function
        If not provided, it guess depending on the data
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
    if plot_func is None:
        order = var.getOrder()
        if var.ndim==1:
            plot_func = stick if quiver else curve
        elif 't' in order:
            plot_func = hov
        elif 'z' in order:
            plot_func = section
        elif order == 'yx':
            plot_func = map
        else:
            plot_func = plot2d

    # Positional arguments
    if plot_func is stick:
        args = vv
    else:
        args = [vv]

    # Default optional arguments
    dict_check_defaults(kwargs,  **DEFAULT_PLOT_KWARGS)

    # Plot
    p = plot_func(*args, **kwargs)

    # Register scalar mappable
    if register_sm:
        sm = p.get_obj('scalar_mappable')
        if sm is None:
            sm = p.axes._gci()
        if sm:
            register_scalar_mappable(p.axes, sm)



def create_map(lon, lat, level=None, axes=None, bathy=None,
               elev=20, azim=-100, add_bathy=True, margin=0, **kwargs):
    """Create a 2D or 3D map with bathy

    Parameters
    ----------
    bathy_<param>
        <param> is passed to :func:`create_bathy_2d` or :func:`create_bathy_3d`
    bathy2d_<param>
        <param> is passed to :func:`create_bathy_2d`
    bathy3d_<param>
        <param> is passed to :func:`create_bathy_3d`
    """
    # Params
    kwbat = kwfilter(kwargs, 'bathy')
    kwbat2d = kwfilter(kwargs, 'bathy2d')
    kwbat3d = kwfilter(kwargs, 'bathy3d')

    # Default optional arguments
    dict_check_defaults(kwargs,  **DEFAULT_PLOT_KWARGS)
    kwargs.update(show=False, close=False)

    # Guess 3D
    if N.isscalar(level) and (level==0 or level=='surf'):
        level = None
    if axes is None and level is not None:
        axes = '3d'
    is3d = (isinstance(axes, Axes3D) or axes=='3d')
    if is3d:
        dict_check_defaults(kwargs,
            drawparallels_linewidth=.1, drawparallels_dashes=[],
            drawmeridians_linewidth=.1, drawmeridians_dashes=[],
            left=0, right=1, bottom=0, top=1, subplot_adjust=True,
        )

    # Map extension
    if isinstance(lon, N.ndarray):
        lon = (lon.min(), lon.max())
    if isinstance(lat, N.ndarray):
        lat = (lat.min(), lat.max())
    if margin:
        box = scalebox(dict(lon=lon, lat=lat), 1+margin)
        lon = box['lon']
        lat = box['lat']

    # Create the map
    m = map(lon=lon, lat=lat, axes=axes, **kwargs)

    # Bathy
    if bathy is not None:
        curv2rect(bathy)
        bathy = bathy(lon=lon, lat=lat)
        if bathy.mean()>0:
            bathy = -bathy
        if level is None:
            level = "bottom"

    # 3D
    if m.is3d:

        # View
        m.axes.set_aspect('auto')
        m.axes.view_init(azim=azim, elev=elev)

        # Set z range
        if level is not None:
            if isinstance(level, (float, int)):
                level = (level, 0)
            elif isaxis(level) or isinstance(level, N.ndarray):
                level = (level[:].min(), 0)
            elif level=='bottom':
                if bathy is not None:
                    level = (bathy.min(), 0)
                else:
                    sonat_warn('Cannot set zlim properly on 3d plot since level=="bottom"'
                        ' and no bathy is provided')
                    level = None
            if level is None:
                level = (-200, 0)
            m.axes.set_zlim(level)

        # Bathy
        if bathy is not None and add_bathy:
            kw = kwbat.copy()
            kw.update(kwbat3d)
            plot_bathy_3d(bathy(lon=lon, lat=lat), m, **kw)

    # 2D bathy
    elif bathy is not None and add_bathy:
        kw = kwbat.copy()
        kw.update(kwbat2d)
        plot_bathy_2d(bathy(lon=lon, lat=lat), m, **kw)

    return m


def plot_bathy_3d(bathy, m, color=land_color, linewidth=0, **kwargs):
    """Plot the bathy as 3D surface

    Parameters
    ----------
    bathy: MV2.array
    m: Map instance
    """
    if bathy.mean()>0:
        bathy = -bathy
    xx = bathy.getLongitude()[:]
    yy = bathy.getLatitude()[:]
    xxb, yyb = N.meshgrid(xx, yy)
    xxb, yyb = m(xxb, yyb)
    bathy = bathy.asma()
    bathy[bathy>0] = 0.
    bathy[bathy<m.axes.get_zlim()[0]] = N.nan
    old_Poly3DCollection = art3d.Poly3DCollection
    art3d.Poly3DCollection = FixedZorderPoly3DCollection
    m.axes.plot_surface(xxb, yyb, bathy, rstride=1, cstride=1,
        color=color, linewidth=linewidth, **kwargs)
    art3d.Poly3DCollection = old_Poly3DCollection

class _DepLabel_(object):
    def __mod__(self, value):
        return '{:.0f} m'.format(abs(value))

def plot_bathy_2d(bathy, m, levels=[-4000, -3000, -2000, -1000, -800, -600,
            -400, -200, -100, -50],
        linewidth=.2, color='.8', linestyle='-', clabel=True,
        clabel_fmt=_DepLabel_(),#'%g m',
        **kwargs):
    """Plot coutours of the bathy on a map
    """
    kwcl = kwfilter(kwargs, 'clabel')
    if bathy.mean()>0:
        bathy = -bathy
    xx = bathy.getLongitude()
    yy = bathy.getLatitude()
    xx, yy = N.meshgrid(xx, yy)
    xx, yy = m(xx, yy)
    cc = m.axes.contour(xx, yy, bathy.asma(), levels=levels, colors=color,
                        linestyles=linestyle, **kwargs)
    if clabel:
        m.axes.clabel(cc, fmt=clabel_fmt, **kwcl)

class FixedZorderPoly3DCollection(Poly3DCollection):
    _zorder = -1
    @property
    def zorder(self):
        return self._zorder
    @zorder.setter
    def zorder(self, value):
        pass


def plot_scattered_locs(lons, lats, depths, slice_type=None, interval=None, plotter=None,
                        lon=None, lat=None, level=None, label='',
                        lon_bounds_margin=.1, lat_bounds_margin=.1,
                        data=None, warn=True, bathy=None, xybathy=None, size=10,
                        linewidth=0.15, color='k', add_profile_line=None, add_bathy=True,
                        fig=None, title=True, register_sm=True,
                        legend=False, **kwargs):
    """Plot scattered localisations

    Parameters
    ----------
    lons: n-D array
    lats: n-D array
    depths: n-D array
    slice_type: one of "3d"/None, "2d", "zonal", "meridional", "horizontal"
        The way to slice the observations.
        "3d"/"2d" are 3D/2D view of all observations.
        Other slices make a selection with a range (``interval``).
    interval: None, tuple of float
        Interval for selecting valid data
        Required if slice_type is not "3d"/None/"2d".
    map_<param>:
        <param> is passed to :func:`create_map`

    Todo
    ----
    Add time support.
    """
    # Inits
    if cdms2.isVariable(data):
        data = data.asma()
    kwmap = kwfilter(kwargs, 'map_')
    kwplt = kwfilter(kwargs, 'plotter_')
    dict_check_defaults(kwmap, **kwplt)
    kwpf = kwfilter(kwargs, 'add_profile_line')
    kwleg = kwfilter(kwargs, 'legend')
    if title is True:
        title = get_long_name(data)

    # Slice type
    if slice_type is None:
        slice_type = "3d"
    else:
        slice_type = str(slice_type).lower()
    valid_slice_types = ['3d', "2d", 'zonal', 'meridional', 'horizonal',
        'bottom', 'surf']
    assert slice_type in valid_slice_types, ('Invalid slice type. '
        'It must be one of: '+', '.join(valid_slice_types))

    # Special depths
    indepths = depths
    if depths=='surf':
        add_profile_line = False
    if (depths=='surf' and slice_type!='surf'): # surface

        depths = N.zeros(len(lons))

    elif (depths=='bottom' and slice_type!='bottom'): # bottom

        if bathy is None:
            raise SONATError('Bathymetry is needed to plot bottom locs with a'
                             ' {} slice type'.format(slice_type))
            depths = grid2xy(bathy, lons, lats)
            if depths.mask.all():
                if warn:
                    sonat_warn('Bathymetry is fully masked at bottom locs. Skipping...')
                return
            if depths.mask.any():
                if warn:
                    sonat_warn('Bathymetry is partly masked at bottom locs. Compressing...')
                valid = depths.mask
                lons = lons[valid]
                lats = lats[valid]
                depths = depths[valid]
                if data is not None:
                    data = data[valid]

    # Pre-slicing
    if (slice_type != '3d' and slice_type!='2d' and
        (slice_type!='surf' or depths!='surf') and
        (slice_type!='bottom' or depths!='bottom')):

        assert interval is not None, ('You must provide a valid '
            '"interval" for slicing scattered locations')
        sliced = slice_scattered_locs(lons, lats, depths, slice_type, interval,
                                      data=data)
        if sliced is None:
            return
        lons, lats, depths, data = sliced

    # Config
    zscattered = (not isinstance(depths, string_types) and not isaxis(depths)
                  and N.ndim(depths)==1 and len(depths)==len(lons))

    # Plotter as Axes
    if isinstance(plotter, P.Axes):
        ax = plotter
        fig = ax.get_figure()
        plotter = None
    elif plotter is None:
        ax = None
    elif isinstance(plotter, Plot):
        ax = plotter.axes
    else:
        raise SONATError('Plotter must be matplotlib Axes instance or '
                         'a vacumm Plot instance')
    if slice_type=='3d':
        if ax is None:
            ax = '3d'
        elif not isinstance(ax, Axes3D):
            sonat_warn("Requesting 3D plot but provided axes are not 3D."
                       " Skipping...")
            axes = None


    # Level for 3d plots
    if slice_type=='3d' and level is None:
        level = depths
        if level.min()==0:
            level = (-200, 0) # Fall back to this interval

    # Get the plotter
    if slice_type in ['3d', "2d", "horizontal", 'surf', 'bottom']: # map

        # Lon/lat
        if lon is None:
            lon = rescale_itv((lons.min(), lons.max()), 1+lon_bounds_margin)
        if lat is None:
            lat = rescale_itv((lats.min(), lats.max()), 1+lat_bounds_margin)

        # Map
        if plotter is None:
            plotter = create_map(lon, lat, level=level, bathy=bathy, add_bathy=add_bathy,
                           fig=fig, axes=ax, **kwmap)
        ax = plotter.axes

        # Projetion
        xx, yy = plotter(lons, lats)

    else:

        xx = lons[:]
        yy = lats[:]

    # Plot params for scatter
    kwargs.update(linewidth=linewidth, s=size)

    # Masking
    if data is not None:
        if data.dtype.char=='?':
            mask = data
            data = None
        else:
            mask = N.ma.getmaskarray(data)
        if zscattered:
            xx = N.ma.masked_where(mask, xx, copy=False)
    elif N.ma.isMA(xx) or N.ma.isMA(yy):
        mask = N.ma.getmaskarray(xx)|N.ma.getmaskarray(yy)
    else:
        mask = None
    if mask is not None and mask.all():
        if warn:
            sonat_warn('All your data are masked')
        return

    # Data kwargs
    if data is not None:
        dict_check_defaults(kwargs, vmin=data.min(), vmax=data.max())

    # 3D
    pp = []
    if slice_type == "3d":

        # Depth labels
        zfmtfunc = lambda x, pos: deplab(x, nosign=True)
        ax.zaxis.set_major_formatter(FuncFormatter(zfmtfunc))

        # Bathy for profile line
        if add_profile_line is None:
            add_profile_line = (indepths != 'surf') and bathy is not None
        if xybathy is not None:
            if add_profile_line and indepths == 'bottom':
                xybathy = indepths
            elif bathy is None:
                if add_profile_line:
                    sonat_warn('Cannot plot profile line without bathymetry')
                add_profile_line = False
            else:
                xybathy = grid2xy(bathy, lons, lats)

        # Scatter plots
        if zscattered: # fully scattered

            # Points
            if data is not None:
                kwargs['c'] = data
            else:
                kwargs['c'] = color
            pp.append(ax.scatter(xx, yy, depths, label=label, **kwargs))

            # Profile lines
            if add_profile_line:
                for ip, (x, y) in enumerate(zip(xx, yy)):
                    plot_profile_line_3d(ax, x, y, xybathy[ip],
                                         zorder=p.get_zorder()-0.01, **kwpf)


        else: # profiles

            for ip, (x, y) in enumerate(zip(xx, yy)):

                isfirst = ip==0

                # Skip fully masked
                if mask is not None:
                    if mask[..., ip].all():
                        if warn:
                            sonat_warn('Profile fully masked')
                        continue

                # Points
                if depths[:].ndim==2:
                    zz = depths[:, ip]
                else:
                    zz = depths
                if mask is not None and mask.ndim==2:
                    zz = N.ma.masked_where(mask[..., ip], zz, copy=False)
                if data is not None:
                    kwargs['c'] = data[..., ip]
                else:
                    kwargs['c'] = color
                pp.append(ax.scatter([x]*len(zz), [y]*len(zz), zz,
                               label=label if isfirst else '', **kwargs))

                # Profile line
                if add_profile_line:
                    plot_profile_line_3d(ax, x, y, xybathy[ip],
                                         zorder=pp[-1].get_zorder()-0.01, **kwpf)

    # 2D
    elif slice_type == "2d":

        # 1D arrays
        if data is not None and data.ndim!=1:
            data = data.mean(axis=0)
        if mask is not None:
            if mask.ndim!=1:
                mask = mask.all(axis=0)
            xx = xx[~mask]
            yy = yy[~mask]
            if data is not None:
                data = data[~mask]

        # Scatter plot
        if data is not None:
            kwargs['c'] = data
        else:
            kwargs['c'] = color
        pp.append(ax.scatter(xx, yy, label=label, **kwargs))


    else:

        raise NotImplementedError('slice_type plot yet not implemented: '+slice_type)

    # Finalise
    if title:
        ax.set_title(title)
    if legend:
        plotter.legend(**kwleg)
    if data is not None and register_sm:
        register_scalar_mappable(ax, pp)

    # TODO: plot of other scattered slice!
    return plotter

def register_scalar_mappable(ax, p):
    if isinstance(p, list):
        for _ in p:
            register_scalar_mappable(ax, _)
    else:
        if not hasattr(ax, '_sonat_scalar_mappables'):
            ax._sonat_scalar_mappables = []
        ax._sonat_scalar_mappables.append(p)

def sync_scalar_mappable_plots_vminmax(ax, symetric=False):
    """Sync min ax max of all the scalar mappable plots for given axes"""
    # Get the axes
    if isinstance(ax, Plot):
        ax = ax.axes
    if not hasattr(ax, '_sonat_scalar_mappables'):
        return

    # Get min and max
    vmin = N.inf
    vmax = -N.inf
    for p in ax._sonat_scalar_mappables:
        vmin = min(p.norm.vmin, vmin)
        vmax = max(p.norm.vmax, vmax)

    # Symetric?
    if symetric:
        vmax = max(abs(vmin), abs(vmax))
        vmin = -vmax

    # Set min and max
    norm = Normalize(vmin, vmax)
    for p in ax._sonat_scalar_mappables:
        p.set_norm(norm)


def plot_profile_line_3d(ax, x, y, bathy, linewidth=.3, linestyle='--', color='.6',
                         size=3, clip=True, **kwargs):
    """Plot a vertical line ending with two point between the bottom and
    the surface"""
    if clip:
        bathy = max(bathy, ax.get_zlim()[0])
    p = ax.plot([x, x], [y, y], [bathy, 0], color=color,
                    linestyle=linestyle, linewidth=linewidth, **kwargs)
    ax.scatter([x, x], [y, y], [bathy, 0], c=color, s=size,
                   zorder=p[0].get_zorder())


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
