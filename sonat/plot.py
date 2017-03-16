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
from cycler import cycler, Cycler
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
from .misc import (slice_gridded_var, vminmax, mask_scattered_locs,
                   rescale_itv, get_long_name)


#: Matplotlib default configuration file
SONAT_DEFAULT_MATPLOTLIBRC =  os.path.join(os.path.dirname(__file__), 'matplotlibrc')

#: Matplotlib user configuration file
SONAT_USER_MATPLOTLIBRC =  'matplotlibrc'

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


#: Default arguments to vacumm plot functions
DEFAULT_PLOT_KWARGS = dict(
    contour_linewidths=.7,
    colorbar_fraction=.1,
    quiver_units='dots',
    quiver_width=1.2,
    fill='contourf',
    contour_zorder=1,
    quiver_norm=3,
    proj='merc',
    cmap='auto',
    drawmeridians_linewidth=rcParams['grid.linewidth'],
    drawmeridians_color=rcParams['grid.color'],
    drawparallels_linewidth=rcParams['grid.linewidth'],
    drawparallels_color=rcParams['grid.color'],
    autoresize='y',
    autoresize_minaspect=.5,
    colorbar_shrink=.8,
#    fillcontinents_zorder=10,
#    figsize=(5, 3), # in matplotlibrc
    show=False,
    close=True,
    )

if rcParams['grid.linestyle'] == '-':
    DEFAULT_PLOT_KWARGS['drawmeridians_dashes'] = ''
    DEFAULT_PLOT_KWARGS['drawparallels_dashes'] = ''
elif rcParams['grid.linestyle'] == ':':
    DEFAULT_PLOT_KWARGS['drawmeridians_dashes'] = [1, 1]
    DEFAULT_PLOT_KWARGS['drawparallels_dashes'] = [1, 1]


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

    return p


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
                        data=None, warn=True, bathy=None, xybathy=None,
                        size=30, color='#2ca02c', linewidth=0.4, edgecolor='k',
                        add_profile_line=None, add_bathy=True,
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
    valid_slice_types = ['3d', "2d", 'zonal', 'merid', 'horiz',
        'bottom', 'surf']
    assert slice_type in valid_slice_types, ('Invalid slice type. '
        'It must be one of: '+', '.join(valid_slice_types))

    # Profiles?
    profiles = (not isinstance(depths, str) and (isaxis(depths) or
                N.shape(depths)!=N.shape(xx) or
                (data is not None and data.ndim==2)))

    # Force some options
    if not profiles or slice_type not in ('3d', 'merid', 'zonal'):
        add_profile_line = False
    elif add_profile_line is None:
        add_profile_line = True

    # Bathymetry
    need_xybathy = int(add_profile_line)
    if depths=='bottom' and slice_type not in ('bottom', '2d'):
        need_xybathy = 2
    if need_xybathy and xybathy is not None:

        if bathy is None and need_xybathy==2: # we really need it
            if warn:
                sonat_warn('Bathymetry is needed at obs locations. Skipping...')
                return
        xybathy = grid2xy(bathy, lons, lats)
        if xybathy.mask.all():
            if warn:
                sonat_warn('Bathymetry is fully masked at bottom locs. Skipping...')
            if need_xybathy==2:
                return

    # Special depths: surf and bottom
    indepths = depths
    if (depths=='surf' and slice_type!='surf'): # surface
        depths = N.zeros(len(lons))
    elif (depths=='bottom' and slice_type not in ('bottom', '2d')): # bottom
        depths = -xybathy
        if interval is not None and N.isscalar(interval[0]):
            interval = (depths+interval[0], depths+interval[1])

    # Numeric coordinates
    xx = lons[:].copy()
    yy = lats[:].copy()
    strdepths = isinstance(depths, str)
    if not strdepths:
        zz = N.array(depths[:], copy=True)

    # Shape
    if data is not None:
        dshape = data.shape
    elif not profiles or strdepths:
        dshape = xx.shape
    elif zz.ndim==2:
        dshape = zz.shape
    else:
        dshape = zz.shape + xx.shape

    # Masking outside interval
    if (slice_type != '3d' and slice_type!='2d' and
        (slice_type!='surf' or depths!='surf') and
        (slice_type!='bottom' or depths!='bottom')):

        assert interval is not None, ('You must provide a valid '
            '"interval" for slicing scattered locations')
        stype = 'horiz' if slice_type in ('surf', 'bottom') else slice_type
        data = mask_scattered_locs(xx, yy, depths, stype, interval,
                                      data=data)
        if data is None:
            return

    # Get the full mask: (np), or (nz, np) for profiles
    # - mask with data
    if data is not None:
        if data.dtype.char=='?':
            mask = data
            data = None
        else:
            mask = N.ma.getmaskarray(data)
    else:
        mask = N.zeros(dshape)
    # - mask with coordinates
    if N.ma.isMA(xx) or N.ma.isMA(yy): # lons/lats
        xymask = N.ma.getmaskarray(xx) | N.ma.getmaskarray(yy)
        mask |= N.resize(mask, dshape)
    if not strdepths and N.ma.isMA(zz): # depths
        if profiles:
            zmask = N.ma.getmaskarray(zz)
            if zz.ndim==1:
                zmask = N.repeat(N.ma.resize(N.ma.getmaskarray(zmask), (-1, 1)),
                           xx.size, axis=1)
            mask |= zmask
        else:
            mask |= N.ma.getmaskarray(zz)
    # - check
    if mask.all():
        if warn:
            sonat_warn('All your data are masked')
        return
    # - mask back
    xymask = mask if mask.ndim==1 else mask.all(axis=0)
    xx = N.ma.masked_where(xymask, xx, copy=False)
    yy = N.ma.masked_where(xymask, yy, copy=False)
    if not strdepths:
        if mask.shape == zz.shape:
            zz = N.ma.masked_where(mask, zz, copy=False)
        elif zz.ndim==1:
            zz = N.ma.masked_where(mask.all(axis=1), zz, copy=False)
    if data is not None:
        data = N.ma.masked_where(mask, data, copy=0)

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


    # Coordinate bounds
    if level is None and slice_type in ['3d', 'zonal', 'merid']:
        if strdepths or zz.min()==0:
            level_min = -200 # Fall back to this min depth
        else:
            level_min = 1.1 * level.min()
        level = (level_min, 0)
    if (lon is None and
        slice_type in ['3d', "2d", "horiz", 'surf', 'bottom', 'zonal']):
        lon = rescale_itv((xx.min(), xx.max()), 1.1)
    if (lat is None and
        slice_type in ['3d', "2d", "horiz", 'surf', 'bottom', 'merid']):
        lat = rescale_itv((yy.min(), yy.max()), 1.1)


    # Get the plotter
    if slice_type in ['3d', "2d", "horiz", 'surf', 'bottom']: # map

        # Map
        if plotter is None:
            plotter = create_map(lon, lat, level=level, bathy=bathy, add_bathy=add_bathy,
                           fig=fig, axes=ax, **kwmap)
        ax = plotter.axes

        # Projection
        xx, yy = plotter(xx, yy)

    else: # sections

        if plotter is None:

            kwplt.update(fig=fig, axes=ax, show=False, close=False)

            if slice_type == 'merid':

                plotter = section(data=None, xaxis=MV2.array(lon, id='lon'),
                                  yaxis=MV2.array(level, id='dep'))
            else:

                plotter = section(data=None, xaxis=MV2.array(lon, id='lat'),
                                  yaxis=MV2.array(level, id='dep'))

        ax = plotter.axes
    axis_bounds = ax.axis()

    # Plot params for scatter
    kwargs.update(linewidth=linewidth, s=size, edgecolor=edgecolor)

    # Data kwargs
    if data is not None:
        dict_check_defaults(kwargs, vmin=data.min(), vmax=data.max())

    # 3D
    pp = []
    if slice_type == "3d":

        # Depth labels
        zfmtfunc = lambda x, pos: deplab(x, nosign=True)
        ax.zaxis.set_major_formatter(FuncFormatter(zfmtfunc))

        # Scatter plots
        if not profiles: # fully scattered

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
                if mask[:, ip].all():
#                    if warn:
#                        sonat_warn('Profile fully masked')
                    continue

                # Points
                if zz.ndim==2:
                    z = depths[:, ip]
                else:
                    z = zz
                z = N.ma.masked_where(mask[:, ip], z, copy=False)
                if data is not None:
                    kwargs['c'] = data[..., ip]
                else:
                    kwargs['c'] = color
                pp.append(ax.scatter([x]*len(z), [y]*len(z), z,
                               label=label if isfirst else '', **kwargs))

                # Profile line
                if add_profile_line:
                    plot_profile_line_3d(ax, x, y, -xybathy[ip],
                                         zorder=pp[-1].get_zorder()-0.01, **kwpf)


    # Horizontal
    elif slice_type in ['2d', 'surf', 'bottom', 'horiz']:

        # Barotropic case
        if slice_type == "2d" and data is not None and data.ndim!=1:
            data = data.mean(axis=0)

        # Scatter plot
        if data is not None:
            kwargs['c'] = data
        else:
            kwargs['c'] = color
        pp.append(ax.scatter(xx, yy, label=label, **kwargs))


    # Sections
    else:

        # X axis data
        if slice_type=='zonal':
            xdata = xx
        else:
            xdata = yy

        # Scatter plots
        if not profiles: # scattered

            if data is not None:
                kwargs['c'] = data
            else:
                kwargs['c'] = color
            pp.append(ax.scatter(xdata, depths, label=label, **kwargs))

        else: # profiles

            for ip, x in enumerate(xdata):

                isfirst = ip==0

                # Skip fully masked
                if mask[:, ip].all():
#                    if warn:
#                        sonat_warn('Profile fully masked')
                    continue

                # Points
                if depths[:].ndim==2:
                    z = zz[:, ip]
                else:
                    z = zz
                z = N.ma.masked_where(mask[:, ip], z, copy=False)
                if data is not None:
                    kwargs['c'] = data[:, ip]
                else:
                    kwargs['c'] = color

                pp.append(ax.scatter([x]*len(z), z,
                               label=label if isfirst else '', **kwargs))

                # Profile line
                if add_profile_line:
                    plot_profile_line_3d(ax, x, -xybathy[ip],
                                         zorder=pp[-1].get_zorder()-0.01, **kwpf)

    # Finalise
    ax.axis(axis_bounds)
    if title:
        ax.set_title(title)
    if legend:
        plotter.legend(**kwleg)
    if data is not None and register_sm:
        register_scalar_mappable(ax, pp)
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
    if bathy > 0:
        bathy = -bathy
    if clip:
        bathy = max(bathy, ax.get_zlim()[0])
    p = ax.plot([x, x], [y, y], [bathy, 0], color=color,
                    linestyle=linestyle, linewidth=linewidth, **kwargs)
    ax.scatter([x, x], [y, y], [bathy, 0], c=color, s=size,
                   zorder=p[0].get_zorder())


def get_color_marker_cycler(colors, markers):
    """Get :class:`cycler.Cycler` to loop over colors and markers"""

    # Individual cyclers
    if colors is None:
        colors = rcParams['axes.prop_cycle']
    elif not isinstance(colors, Cycler):
        if not isinstance(colors, list):
            colors = [colors]
        colors = cycler(color=colors)
    if not isinstance(markers, Cycler):
        if not isinstance(markers, list):
            markers = [markers]
        markers = cycler(marker=markers)

    # Sync them
    if len(colors) < len(markers):
        colors = colors * (len(markers) / len(colors) + 1)
        colors = colors[:len(markers)]
    elif len(colors) > len(markers):
        markers = markers * (len(colors) / len(markers) + 1)
        markers = markers[:len(colors)]

    # Combine them
    cyc = colors + markers

    # Return generator
    return cyc()


