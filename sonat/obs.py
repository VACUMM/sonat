"""
Observations
"""
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

from collections import OrderedDict
import os
import re
from six import string_types
import cdms2
import MV2
import numpy as N
from cycler import cycler
from matplotlib.markers import MarkerStyle
from vcmq import (grid2xy, regrid2d, ncget_lon, ncget_lat,
    ncget_time, ncget_level, ArgList, ncfind_obj, itv_intersect, intersect,
    MV2_axisConcatenate, transect, create_axis, regrid1d, isaxis,
    dicttree_get, dict_check_defaults, meshgrid, P, kwfilter)

from .__init__ import sonat_warn, SONATError, BOTTOM_VARNAMES, get_logger
from .misc import (xycompress, _Base_, _XYT_, check_variables, _NamedVariables_,
                   rescale_itv, get_long_name, split_varname)
from .pack import default_missing_value
from .stack import Stacker, _StackerMapIO_
from .plot import plot_scattered_locs, sync_scalar_mappable_plots_vminmax

npy = N

RE_PERTDIR_MATCH = re.compile(r'[+\-][xy]$').match

OBS_PLATFORM_TYPES = {}

def register_obs_platform(cls, warn=True, replace=True):
    """Register a new observation platform type"""

    # Check class
    if not issubclass(cls, ObsPlatformBase):
        raise SONATError('Platform cannot be registered: '
            'it must be a subclass of ObsPlatformBase')

    # Check platform_type existence and type
    if (not hasattr(cls, 'platform_type') or
            not isinstance(cls.platform_type, basestring)):
        raise SONATError('Platform cannot be registered: no attribute '
            '"platform_type"')

    # Check platform_type value
    if cls.platform_type in OBS_PLATFORM_TYPES:
        if warn:
            msg = 'Plaform already registered: {}'.format(cls.platform_type)
            if replace:
                msg = msg + '. Replacing it...'
            sonat_warn(msg)
        if not replace:
            return

    OBS_PLATFORM_TYPES[cls.platform_type] = cls

def get_obs_platform(platform_type, *args, **kwargs):
    """Get the class or instance for a given platform type

    If extra arguments are passed, they are used to create a instance
    which is then return. Otherwise, the class is return.
    """

    # Get the class
    if platform_type not in OBS_PLATFORM_TYPES:
        raise SONATError('Observation platform type not registered: '+platform_type)
    cls = OBS_PLATFORM_TYPES[platform_type]

    # Instance
    if args or kwargs:
        return cls(*args, **kwargs)
    return cls

def load_obs_platform(platform_type, pfile, varnames=None, name=None, **kwargs):
    """Load an aobservation platform instance from its type and its file"""
    obs = get_obs_platform(platform_type, pfile, varnames=varnames, **kwargs)
    obs.name = name
    return obs


class _ObsBase_(_XYT_):

    def set_cached_plot(self, var_name, slice_type, slice_loc, obj):
        """Cache a plotting object knowing the slice and the variable name"""
        obj._sonat_plot = var_name, slice_type, slice_loc
        self.plot_cache.setdefault(var_name, {}).setdefault(
            slice_type, {}).setdefault(slice_loc, obj)

    def get_cached_plot(self, var_name, slice_type, slice_loc, default=None):
        """Get a cached plot object knowing the slice and the variable name

        Return
        ------
        Map, axes, None
            None is returned if not cached
        """
        return self.plot_cache.get(var_name, {}).get(slice_type, {}
                    ).get(slice_loc, default)

    def set_plot_cache(self, cache):
        assert isinstance(cache, dict), 'Plot cache must a dictionary'
        self._plot_cache = cache

    def get_plot_cache(self):
        if not hasattr(self, '_plot_cache'):
            self._plot_cache = {}
        return self._plot_cache

    def del_plot_cache(self):
        if hasattr(self, '_plot_cache'):
            del self._plot_cache

    plot_cache = property(fget=get_plot_cache, fset=set_plot_cache,
                          fdel=del_plot_cache, doc='Plot cache')

    def get_cached_plots(self):
        """Get the list of all cached plots"""
        plotters = []
        for pv in self.plot_cache.values():
            for pt in pv.values():
                plotters.extend(pt.values())
        return plotters

    def save_cached_plot(self, plotter, figpat, **subst):
        """Save the figure of the given plotter

        Parameters
        ----------
        plotter: vacumm or mpl plot instance
        figpat: string
            Figure file pattern that support substitutions
        \**subst: dict
            Dictionary of string to substitute in figpat

        Return
        ------
        figs
            dict(var_name={slice_type: {slice_loc:figfile}})
        """
        if not hasattr(plotter, '_sonat_plot'):
            raise SONATError('This plot has not been cached')
        var_name, slice_type, slice_loc = plotter._sonat_plot
        subst = subst.copy()
        subst.update(**locals())
        figfile = figpat.format(**subst)
        plotter.savefig(figfile)
        self.created(figfile)
        return dict(var_name={slice_type: {slice_loc:figfile}})


    def save_cached_plots(self, figpat, **subst):
        """Save all cached plots"""
        figs = {}
        for plotter in self.get_cached_plots():
            figs.update(self.save_cached_plot(plotter, figpat, **subst))
        return figs


    def mark_cached_plot_legend(self, plotter):
        """Mark a plot to indicate it has a pending legen to plot
        with :meth:`add_cached_plots_legend`
        """
        plotter._sonat_legend = True

    def add_cached_plots_legend(self, **kwargs):
        for plotter in self.get_cached_plots():
            if hasattr(plotter, '_sonat_legend'):
                plotter.legend(**kwargs)


    def get_level(self, bathy2d=None, margin=0, zmax=0):
        """Get the numeric depth at observation locations

        .. warning:: If observations are at the bottom, the method requires
            the ``bathy2d`` argument to estimate the bottom depths using
            bilinear interpolation.
        """
        depths = self.get_num_depths(bathy2d)
        level = rescale_itv((depths[:].min(), depths[:].max()), factor=margin+1)
        if level[1]>zmax:
            level = level[:1] + (zmax,) + level[2:]
        return level


class ObsPlatformBase(Stacker, _ObsBase_, _NamedVariables_, ):

    name = None

    @property
    def platform_type(self):
        raise SONATError('"platform_type" must declared as a string attribute')

    def __init__(self, *args, **kwargs):
        raise SONATError('This method must be overwritten')

    def project_model(self, *args, **kwargs):
        raise SONATError('This method must be overwritten')

class NcObsPlatform(ObsPlatformBase):
    """Generic observation platform class


    Attributes
    ----------
    pshape: string
        Plateform shape
    varnames: strings
        Variable names WITHOUT the error prefix
    """

    nc_error_suffix = '_error'
    nc_mobility_name = 'mobility'
    platform_type = 'generic'

    def __init__(self, ncfile, varnames=None, logger=None, norms=None,
            time=None, lon=None, lat=None, level=None, name=None,
            singlevar=False, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Init parameters
        self.ncfile = ncfile
        self.varnames = ArgList(varnames).get() if varnames else varnames
        self.time = time
        self.lon = lon
        self.lat = lat
        self.level = level
        self.errors = OrderedDict()
        self._orig = {}
        self.name = name

        # Load positions and variables
        self.load(singlevar=singlevar)

        # Init stacker
        Stacker.__init__(self, self.errors.values(), logger=False, means=False,
            norms=None if isinstance(norms, dict) else norms)

        # Named norms
        if norms and isinstance(norms, dict):
            self.set_named_norms(norms)


    def load(self, singlevar=False, **kwargs):
        """Read mobility, errors, depth and time in open netcdf file"""
        # Open
        f = cdms2.open(self.ncfile)

        # List of error variables
        if self.varnames is None:

            # List from file without error suffix
            self.varnames = [vname.rstrip(self.nc_error_suffix)
                for vname in f.listvariables()
                    if vname.endswith(self.nc_error_suffix)]
        else:

            # Remove error suffix
            self.varnames = [vname.rstrip(self.nc_error_suffix)
                for vname in self.varnames]

            # Keep valid names
            self.varnames = [vname for vname in self.varnames
                if vname+self.nc_error_suffix in f.listvariables()]

        if not self.varnames:
            raise SONATError(('No valid error variable with suffix "{0.nc_error_suffix}"'
                ' in file: {0.ncfile}').format(self))
        if singlevar:
            self.varnames = self.varnames[:1]

        # Reference variable
        fs = f[self.varnames[0] + self.nc_error_suffix]

        # Inspect
        grid = fs.getGrid()
        order = fs.getOrder()
        kwread = {}
        if grid and len(grid.shape)==2: # Structured grid

            self.pshape = 'gridded'
            kwread = dict(time=self.time, lat=self.lat, lon=self.lon, level=self.level)
            if '-' in order:
                raise SONATError("There are unkown dimensions in your gridded variable. "
                    "Current order: "+order)
            order = order[:-2]


        else: # Unstructured grid

            mask = N.ma.nomask #zeros(fs.shape[-1], '?')
            self.pshape = 'xy'
            paxis = fs.getAxis(-1)
            kwread = {}

            # Check order
            order = order[:-1]
            if '-' in order:
                raise SONATError("There are unkown dimensions in your gridded variable. "
                    "Current order: "+order)

            # Lon/lat selection
            if grid:
                kwread.update(lat=self.lat, lon=self.lon)
            else:
                lons = ncget_lon(f)
                if lons is None:
                    raise SONATError("Longitudes not found")
                if lons.getAxis(0).id is not paxis.id:
                    raise SONATError("Longitudes dimension must be the same"
                                     " as last variable dimension")
                lats = ncget_lat(f)
                if lats is None:
                    raise SONATError("Latitudes not found")
                if lats.getAxis(0).id is not paxis.id:
                    raise SONATError("Latitudes dimension must be the same"
                                     " as last variable dimension")
                lons = lons.asma()
                lats = lats.asma()
                if self.lon:
                    mask |= lons < self.lon[0]
                    mask |= lons > self.lon[1]
                if self.lat:
                    mask |= lats < self.lat[0]
                    mask |= lats > self.lat[1]

            # Vertical dimension
            if 'z' in order: # 1D axis
                kwread.update(level=self.level)
            else: # Aux axis or variable
                self.depths = ncget_level(f)
                if self.depths is not None:
                    if (self.depths[:].ndim==1 and
                            self.depths.getAxis(0).id is paxis.id):
                        raise SONATError("Scattered Z dimension not yet supported")
                        self.pshape = self.pshape + 'z'
                    depths = self.depths.asma()
                    if self.level:
                        mask |= depths < self.level[0]
                        mask |= depths > self.level[1]

            # Time selection
            if 't' in order: # 1D axis
                kwread.update(time=self.time)
            else: # Aux axis
                self.times = ncget_time(f)
                if self.times is not None:
                    if self.time:
                        times = create_time(times[:], times.units)[:]
                        mask |= times < reltime(self.time[0], times.units)
                        mask |= times > reltime(self.time[1], times.units)
                    self.pshape = self.pshape + 't'

            # Check mask
            if N.ma.isMA(mask):
                mask = mask.filled(False)
            if mask.all():
                SONATError("All your observation data are masked")

        # Read
        # - errors
        for vname in self.varnames:
            self.errors[vname] = f(vname + self.nc_error_suffix, **kwread)
        self._sample = sample = self.errors[vname]
        order = sample.getOrder()
        self.axes = []
        # - lon/lat
        if grid:
            self.lons = sample.getLongitude()[:]
            self.lats = sample.getLatitude()[:]
        else:
            self.lons = lons
            self.lats = lats
        # - times
        if 't' in order:
            self.times = sample.getTime()
            if self.times is not paxis:
                self.axes.append(self.times)
        else:
            self.times = None
        # - depths
        if 'z' in order:
            self.depths = sample.getLevel()
            self.axes.append(self.depths)
        elif (hasattr(self, 'depths') and self.depths is not None
              and 'z' not in self.pshape): # depths that vary with position
            self.axes.append(self.depths)
        elif hasattr(f, 'depth'): # depth from file attribute
            self.depths = f.depth
            if not isinstance(self.depths, basestring):
                raise SONATError("Depth must be a string if specified as an attribute")
        elif fs.id in BOTTOM_VARNAMES: # 2D: bottom
            self.depths = 'bottom'
        else: # surf by default
            self.depths = 'surf'
        # - mobility
        if self.nc_mobility_name in f.listvariables():
            self.mobility = f(self.nc_mobility_name, **kwread)
        elif hasattr(f, self.nc_mobility_name): # as an attribute
            self.mobility = getattr(f, self.nc_mobility_name)
        elif self.pshape=='gridded': # don't move grid for the moment
            self.mobility = 0.
        else:
            sonat_warn('No mobility variable found in obs file. '
                'Setting it to 1 everywhere.')
            shape = grid.shape if grid else lons.shape
            self.mobility = MV2.ones(shape, 'i')
        self.mobile = (N.ma.asarray(self.mobility)==1).all()

        # Compression to remove masked points
        if self.pshape != 'gridded' and mask.any():
            valid = ~mask
            self.lons = xycompress(valid, self.lons)
            self.lats = xycompress(valid, self.lats)
            for vname, var in self.errors.items():
                self.errors[vname] = xycompress(valid, var)
            if self.times is not None and self.times not in self.axes:
                self.times = xycompress(valid, self.times)
            if (self.depths is not None and not isaxis(self.depths) and
                    not isinstance(self.depths, basestring)):
                self.depths = xycompress(valid, self.depths)

        # Name
        if self.name is None:
            if hasattr(f, 'platform_name'): # name from attribute
                self.name = f.platform_name
            else: # name from file name
                self.name = os.path.splitext(os.path.basename(self.ncfile))[0]
        f.close()

    @property
    def platform_name(self):
        return self.name

    def get_suffixed_varname(self, vname):
        """Get the name of var optionally with a depth suffix"""
        if hasattr(self, 'id'):
            vname = vname.id
        if isinstance(self.depths, basestring):
            vname = vname + '_' + self.depths
        return vname

    @property
    def suffixed_varnames(self):
        return [self.get_suffixed_varname(vname) for vname in self.varnames]

    @property
    def ndim(self):
        return self._sample.ndim

    @property
    def nsdim(self):
        return 1 + int(self.pshape=='gridded')

    @property
    def shape(self):
        return self._sample.shape

    @property
    def grid(self):
        return self._sample.getGrid()

    @property
    def is_surf(self):
        return self.depths is 'surf'

    @property
    def is_bottom(self):
        return self.depths is 'bottom'

    @property
    def is_zscattered(self):
        return 'z' in self.pshape

    @property
    def is_tscattered(self):
        return 't' in self.pshape

    @property
    def is_gridded(self):
        return self.pshape is 'gridded'

    @property
    def mask(self):
        return self.get_mask()

    def get_mask(self, scattered=False):
        """Get the mask of data

        Parameters
        ----------
        scattered: bool
            Flatten in X and Y, for gridded platforms only.
        """
        mask = N.ma.getmaskarray(self._sample)
        if scattered and self.is_gridded:
            mask = mask.reshape(mask.shape[:-2] + (-1,))
        return mask

    @property
    def lons1d(self):
        """Longitudes raveled in space"""
        if not self.is_gridded:
            return self.lons
        self._check_lonslats1d_()
        return self._lons1d

    @property
    def lats1d(self):
        """Latitudes raveled in space"""
        if not self.is_gridded:
            return self.lats
        self._check_lonslats1d_()
        return self._lats1d

    def _check_lonslats1d_(self):
        if not hasattr(self, '_lons1d'):
            _lons1d, _lats1d = meshgrid(self.lons, self.lats)
            self._lons1d = _lons1d.ravel()
            self._lats1d = _lats1d.ravel()

    def xy_ravel_var(self, var):
        """Ravel variable in space"""
        if var is None:
            return
        if not self.is_gridded:
            return var
        varm = var.reshape(var.shape[:-2] + (-1, ))
        if not cdms2.isVariable(var):
            return varm
        varo = MV2.array(varm, copy=0, attributes=varm.attributes, id=var.id)
        for i, ax in enumerate(var.getAxisList()[:-2]):
            varo.setAxis(i, ax)
        return varo


    def check_variables(self, searchmode='ns'):
        """Check that all input variables have known properties

        See also
        --------
        :func:`sonat.misc.check_variables`
        """
        check_variables([pack.input for pack in self], searchmode=searchmode)

    def get_valid_varnames(self, varnames=None):
        """Filter out var names that are not valid or not equal to 'locations'

        If names are explicit numeric arrays, they are unchanged
        """
        if varnames is None:
            return None
        if isinstance(varnames, string_types):
            varnames = [varnames]
        for vn in list(varnames):
            if vn is True:
                for svn in self.varnames:
                    if svn not in varnames:
                        varnames.extend(svn)
        return filter(lambda x: x in self.varnames or x=='locations', varnames)

    def get_error(self, vname, scattered=False):
        """Get an error variable

        Parameters
        ----------
        vname: string
            Name of the variable
        scattered: bool
            Flatten in X and Y, for gridded platforms only.
        """
        var = self.errors[vname]
        if scattered and self.is_gridded:
            var = var.reshape(var.shape[:-2] + (-1,))
        return var

    def set_named_norms(self, *anorms, **knorms):
        """Set norms by variable names

        Note
        ----
        The :class:`~sonat.pack.Packer` instance of variables that were
        not normed are returned.

        Example
        -------
        >>> obs.set_named_norms(dict(temp=2.), sal=.8)
        """
        dnorms = dict(*anorms, **knorms)
        notnormed = []
        restacked = False
        for packer in self:
            for varname, norm in dnorms.items():
                if packer.id and packer.id.split('_')[0] == varname:
                    packer.norm = norm
                    restack = True
                    break
            else:
                notnormed.append(packer)
        if restack:
            self._core_stack_()
        return notnormed


    def get_named_minmax(self):
        """Return (min,max) for each variables in a dict"""
        mm = {}
        for i, name in enumerate(self.varnames):
            if name:
                mm[name] = self[i].data.min(), self[i].data.max()
        return mm


    def activate_xy_pert(self, pert, zonal, index):
        """

        Parameters
        ----------
        pert: float
            Positive or nagative position perturbation in meters
        zonal: bool
            Zonal perturbation? Else meridional.
        index: int
            Index of a mobile observation
        """
        if not self.mobile:
            sonat_warn('This platform cannot be moved')
            return

        if direction:

            # Parse direction
            direction = str(direction).lower()
            if not RE_PERTDIR_MATCH(direction):
                raise SONATError("The direction of perturbation argument must"
                    " be of the form: '{+|-}{x|y}'")
            isx = direction[1]=='x'
            sign = 1 if direction[0]=='-' else -1

            # X and Y perturbation values
            if zonal:
                lat = self.lats.mean()
                pert = m2deg(pert, lat)
                cname = 'lons'
            else:
                cname = 'lats'
            coords = getattr(self, cname) # array of coordinates

            # Reset
            self.deactivate_pert()

            # Save original coordinates values
            self._orig[cname] = coords.copy()

            # Change coordinates in place
            coords[index] += pert

        elif self._orig: # back to orig

            self.deactivate_pert()


    def deactivate_xy_pert(self):
        if not self.mobile:
            return
        for cname in 'lons', 'lats':
            if cname in self._orig:
                getattr(sel, cname)[:] = self._orig[cname]
                del self._orig[cname]

    def get_xy_pert_indices_iter(self):
        if not self.mobile:
            return iter([])
        if not hasattr(self, '_xy_pert_indices'):
            self._xy_pert_indices = N.where(self.mobility>=1)
        return iter(self._xy_pert_indices)

    def init_xy_pert_arrays(self):
        """Get the results for a spatial perturbation sensitivity test

        Return
        ------
        array
            First order zonal derivative
        array
            First order meridional derivative
        array
            Second order zonal derivative
        array
            Second order meridional derivative
        """
        if not self.mobile:
            sonat_warn('This platform cannot be moved')
            return

        mobile = self.mobility>=1
        xder1 = xycompress(mobile, self._sample, id='xderiv1',
            long_name='First order zonal derivative of variance')
        yder1 = xder1.clone()
        yder1 = 'xderiv1'
        yder1.long_name = 'First order meridional derivative of variance'
        xder2 = xder1.clone()
        xder2.id = 'xderiv2'
        xder2.long_name = 'Second order zonal derivative of variance'
        yder2 = xder2.clone()
        yder2.id = 'yderiv2'
        yder2.long_name = 'Second order meridional derivative of variance'

        return xder1, yder1, xder2, yder2

    def project_model(self, var, checkid=True):
        """Project model variables to observations positions"""
        # List of variables
        al = ArgList(var)

        # Loop on variables
        out = []
        for i, var in enumerate(al.get()):

            # Check id
            if checkid and var.id not in self.suffixed_varnames:
                sonat_warn('Variable id "{}" not found on this observation platform'.format(
                    var.id))

            # Order
            order = var.getOrder()

            # Gridded and scattered
            if self.pshape=='gridded':
                var = regrid2d(var, self.grid, method='bilinear')
            else:
                kw = {}
                if 't' in self.pshape:
                    if 't' not in order:
                        raise SONATError("Model variables must have a time axis")
                    kw['times'] = self.times
                if 'z' in self.pshape:
                    if 'z' not in order:
                        raise SONATError("Model variables must have a depth axis")
                    kw['depths'] = self.depths

                if i==0:
                    outaxis = create_axis(self.lons.shape, id='point',
                        long_name='Observation point')

                var = transect(var, lons=self.lons, lats=self.lats,
                    outaxis=outaxis, **kw)


            # Auxilary axes
            for axis in self.axes:
                kw = {}
                if axis[:].ndim==2: # varying depths
                    kw.update(axis=-2, iaxi=0)
                var = regrid1d(var, axis, method='linear', **kw)

            out.append(var)

        return al.put(out)

    def plot(self, variables=None,
             full3d=True, full2d=True, surf=None, bottom=None, horiz_sections=None,
             zonal_sections=None, merid_sections=None,
             lon=None, lat=None, level=None,
             lon_interval_width=0.1, lat_interval_width=0.1, dep_interval_width=0.1,
             bathy=None, m=None, fig=None, close=True,
             savefig=True, add_bathy=True, add_profile_line=None,
             vmin=None, vmax=None, cmap=None,
             figpat='sonat.obs.{platform_type}_{platform_name}_{var_name}_{slice_type}_{slice_loc}.png',
             reset_cache=True, label=None, legend=True, title=True,
             zorder=2.5, sync_vminmax=True, subst={},
             **kwargs):
        """Plot observations locations or data

        Parameters
        ----------
        variables: strings, arrays
            May be variable names, or array that are comptible
            with this platform.
            The special value "locations" just plot the locations.
        lon_interval_width: float in degrees
            Longitude interval width to collect obs points
        lat_interval_width: float in degrees
            Latitude interval width to collect obs points
        dep_interval_width: float in meters
            Depth interval width to collect obs points
        """
        # Inits
        self.verbose('Plotting observations for '+self.name)
        figs = OrderedDict()
        platform_name = self.platform_name
        platform_type = self.platform_type
        if reset_cache:
            del self.plot_cache
        if label is None:
            label = self.platform_name
        kwleg = kwfilter(kwargs, 'legend')

        # Data
        if variables is None:
            variables = ['locations']
        varspecs = []
        if not isinstance(variables, list):
            variables = [variables]
        for var in variables:
            if not isinstance(var, string_types):
                var_name = var.id
            elif var == 'locations':
                var_name = var
                var = self.get_mask(scattered=True)
            else:
                var_name = var
                var = self.get_error(var, scattered=True)
            varspecs.append((var_name, var))

        # Bathy for profile line
        xybathy = None
        if add_profile_line is None:
            add_profile_line = not self.is_zscattered
        if not (full3d or zonal_sections or merid_sections):
            add_profile_line = False
        if add_profile_line:
            if bathy is None:
                sonat_warn('Bathymetry is needed to plot profile lines')
                add_profile_line = False
            else:
                xybathy = self.get_bathy(bathy)

        # Bounds
        if (level is None and (self.depths!='bottom' or bathy is not None) and
            (full3d or zonal_sections or merid_sections)):
            level = self.get_level(bathy)[0], 0

        # 2 and 3D full plots
        slice_type = 'map'
        slice_locs = []
        if full3d:
            slice_locs.append('3d')
        if full2d:
            slice_locs.append('2d')
        for slice_loc in slice_locs:
            self.debug(' Slice: {} / {}'.format(slice_type, slice_loc))

            for var_name, var in varspecs:
                var_name = split_varname(var_name)[0] # no suffixes
                varname = var_name
                self.debug('  Var: ' + var_name)
                long_name = get_long_name(var, var_name)

                # Local args
                default_plotter = dicttree_get(m, var_name, slice_type, slice_loc)
                this_plotter = self.get_cached_plot(var_name, slice_type, slice_loc,
                    default=default_plotter)
                this_fig = dicttree_get(fig, var_name, slice_type, slice_loc)
                if this_fig is None and this_plotter is None: # no active plot
                    this_fig = 'new'
                this_add_bathy = dicttree_get(add_bathy, var_name, slice_type, slice_loc)
                this_vmin = dicttree_get(vmin, var_name, slice_type, slice_loc)
                this_vmax = dicttree_get(vmax, var_name, slice_type, slice_loc)
                this_cmap = dicttree_get(cmap, var_name, slice_type, slice_loc)
                this_title = dicttree_get(title, var_name, slice_type, slice_loc)
                if this_title is True:
                    this_title = long_name

                kw = kwargs.copy()
                if full3d and slice_loc=="3d":
                    this_add_profile_line = dicttree_get(add_profile_line, var_name,
                        slice_type, slice_loc)
                    kw.update(level=level, add_profile_line=this_add_profile_line,
                              xybathy=xybathy)

                # Generic scatter plot
                this_plotter = plot_scattered_locs(
                                    self.lons1d, self.lats1d, self.depths,
                                    slice_type=slice_loc,
                                    data=self.xy_ravel_var(var),
                                    plotter=this_plotter, fig=this_fig,
                                    warn=False, bathy=bathy, label=label,
                                    vmin=this_vmin, vmax=this_vmax, cmap=this_cmap,
                                    lon=lon, lat=lat,
                                    legend=False, title=this_title,
                                    add_bathy=this_add_bathy, zorder=zorder,
                                    **kw)
                if this_plotter is None:
                    continue

                # Sync vmin/vmax with other plots
                if sync_vminmax:
                    sync_scalar_mappable_plots_vminmax(this_plotter, sync_vminmax)

                # Cache
                self.set_cached_plot(var_name, slice_type, slice_loc, this_plotter)

                # Legend
                if legend:
                    self.mark_cached_plot_legend(this_plotter, **kwleg)
                    if legend is True:
                        plotter.legend(**kwleg)

                # Save
                if savefig:
                    figs.update(self.save_cached_plot(this_plotter, figpat, **subst))



        # TODO: plot other slices in obs
        return figs


    def get_bathy(self, bathy2d=None):
        """Get the bathymetry interpolated to XY locations"""
        if hasattr(self, '_bathy'):
            return self._bathy
        if bathy2d is None:
            return None
        self._bathy = grid2xy(bathy2d, self.lons1d, self.lats1d)
        return self._bathy

    def get_num_depths(self, bathy2d=None):
        """Get the depths as numeric array"""
        if not isinstance(self.depths, string_types):
            return self.depths
        if self.depths == 'surf':
            return self.lons*0.
        if bathy2d is None:
            self.error("Can't get numeric bottom depths without bathymetry")
        return self.get_bathy(bathy2d)

class ObsManager(_Base_, _StackerMapIO_, _ObsBase_):
    """Class to manage several observation platform instances"""

    def __init__(self, input, logger=None, norms=None, syncnorms=True,
            missing_value=default_missing_value, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Load platforms
        obsplats = _StackerMapIO_._load_input_(self, input)
        self.obsplats = []
        self.inputs = []
        for obsplat in obsplats:
            if not isinstance(obsplat, NcObsPlatform):
                raise SONATError('ObsManager must be initialised with a single or '
                    'list of NcObsPlatform instances')
            obsplat.set_missing_value(missing_value)
            self.obsplats.append(obsplat)
            self.inputs.append(obsplat.inputs)
        self._missing_value = missing_value

        # Named norms
        if norms and isinstance(norms, dict):
            self.set_named_norms(norms)

        # Synchronise norms
        self._norm_synced = False
        if syncnorms:
            self.sync_norms()

        # Stack stacked data
        self._core_stack_()
        self.splits = npy.cumsum([obs.stacked_data.shape[0]
            for obs in self.obsplats[:-1]])

    def _core_stack_(self):
        self.stacked_data = npy.asfortranarray(
            npy.concatenate([obs.stacked_data for obs in self.obsplats], axis=0))

    def __len__(self):
        return len(self.obsplats)

    def __getitem__(self, key):
        return self.obsplats[key]

    def __iter__(self):
        for obs in self.obsplats:
            yield obs

    def set_missing_value(self, missing_value):
        if missing_value != self._missing_value:
            for obsplats in self.obsplats:
                obsplats.set_missing_value(missing_value)
            self._missing_value = missing_value

    def get_missing_value(self):
        return self._missing_value

    missing_value = property(fget=get_missing_value, fset=set_missing_value)
    fill_value = missing_value

    @property
    def varnames(self):
        vv = []
        for obs in self:
            vv.extend(obs.varnames)
        return list(set(vv))

    @property
    def suffixed_varnames(self):
        vv = []
        for obs in self:
            vv.extend(obs.suffixed_varnames)
        return list(set(vv))

    def check_variables(self, searchmode='ns'):
        """Check that all input variables have known prorties

        See also
        --------
        :func:`sonat.misc.check_variables`
        """
        for obs in self:
            obs.check_variables(searchmode)

    @property
    def lons(self):
        if not hasattr(self, '_lons'):
            self._lons = N.ma.concatenate([obs.lons[:].ravel() for obs in self
                if obs.lons is not None])
        return self._lons

    @property
    def lats(self):
        if not hasattr(self, '_lats'):
            self._lats = N.ma.concatenate([obs.lats[:].ravel() for obs in self
                if obs.lats is not None])
        return self._lats

    @property
    def times(self):
        if not hasattr(self, '_times'):
            times = [obs.times for obs in self if obs.times is not None]
            if not times:
                self._times = None
            else:
                self._times = MV2_axisConcatenate(times)
        return self._times

    @property
    def depths(self):
        """Dict of depths for each variable across all platforms"""
        if not hasattr(self, '_depths'):

            # Loop on platforms
            self._depths = {}
            for obs in self:
                if isinstance(obs.depths, basestring): # 2D
                    depths = obs.depths
                elif obs.depths is not None: # 3D
                    depths = '3d'
                else:
                    depths = None
                for varname in obs.varnames:
                    vdepths = self._depths.setdefault(varname, [])
                    if depths is not None and depths not in vdepths:
                        vdepths.append(depths)

            # Finalisation for each variable
            for varname in self.varnames:
                if not self._depths[varname]:
                    self._depths[varname] = None
                else:
                    self._depths[varname] = tuple(self._depths[varname])

        return self._depths


    def get_num_depths(self, bathy2d=None):
        """Get the depths as a numeric flat array"""
        return N.ma.concatenate([obs.get_num_depths(bathy2d) for obs in self])

    @property
    def has_bottom(self):
        return any([obs.is_bottom for obs in self])


#    def _get_R(self):
#        if not hasattr(self, '_R'):
#            self._R = N.asfortranarray(N.diag(self.stacked_data))
#        return self._R


    def set_norms(self, norms):
        """Change the norms of all abservation plateforms

        Parameters
        ----------
        norms: list of list

        Example
        -------
        >>> obsmanager.set_norms([[2., 4.], [1., .4, 6]])
        """
        norms = self.remap(norms, reshape=True)
        for obs, norm in zip(self, norms):
            obs[i].set_norms(norm)
        self._core_stack_()

    def set_named_norms(self, *anorms, **knorms):
        """Set norms by variable names

        Note
        ----
        The :class:`~sonat.pack.Packer` instance of  allvariables that were
        not normed is returned.

        Example
        -------
        >>> obsmanager.set_named_norms(dict(temp=2.), sal=.8)
        """
        dnorms = dict(*anorms, **knorms)
        notnormed = []
        for obs in self:
            notnormed.extend(obs.set_named_norms(dnorms))
        self._core_stack_()
        return notnormed

    def sync_norms(self, force=True):
        """Synchronise norms accross platforms between variables with the same name

        Parameters
        ----------
        force: bool
            Force sync even if it seems at has already been synced.
        """
        if not force and self._norm_synced:
            return False

        # Renorm for each variable
        for varname in self.varnames:

            # Unified norm
            psize = 0.
            sqnorms = 0.
            packers = []
            for obs in self.obsplats:
                if varname in obs.varnames:
                    i = obs.varnames.index(varname)
                    packers.append(obs[i])
                    psize += obs[i].psize
                    sqnorms += obs[i].psize * obs[i].norm**2
            if len(packers)==1: # nothing to unify
                continue
            norm = N.sqrt(sqnorms / psize)

            # Set it (same as self.set_named_norms({'varname':norm})
            for packer in packers:
                packer.set_norm(norm)

        self._norm_synced = True
        if hasattr(self, 'stacked_data'):
            self._core_stack_()
        return True

    def sync_scatter_plots_vminmax(self):
        """Sync min and max of all scatter plots that are on the same axes"""
        for ax in self.get_cached_plots():
            sync_scatter_plots_vminmax(ax)


    def get_model_specs(self):
        """Get specifications for reading model

        Return
        ------
        dict
            Keys:

            - lon: longitude interval
            - lat: latitude interval
            - depths: depths as a dict for each variable within 'surf', 'bottom', '3d'
            - varnames: variable names
        """
        specs = self.get_seldict()
        specs['depths'] = self.depths
        specs['varnames'] = self.varnames
        return specs


    def restack(self, input, scale='norm'):

        # Check input data
        inputs = self.remap(input, reshape=False)

        # Stack stacks
        stacks = [self[istack].restack(unstacked, scale=scale)
            for istack, unstacked in enumerate(inputs)]

        # Check record length (first axis)
        nr1 = stacks[0].size / self[0].ns
        if len(inputs)>1:
            for i, d in enumerate(stacks[1:]):
                i += 1
                nr2 = d.size / self[i].ns
                if stacks[0].ndim != stacks[i].ndim:
                    if stacks[0].ndim==1:
                        SONATError('Variable {} must not have a record dimension '.format(i))
                    else:
                        SONATError('Variable {} must have a record dimension '.format(i))
                elif nr1!=nr2:
                    raise SONATError(('Record length of variable {i} ({nr2}) '
                        'different from that of first variable ({nr1})').format(**locals()))

        # Stack
        stacked_data = npy.asfortranarray(npy.concatenate(stacks, axis=0))

        return stacked_data

    def unstack(self, sdata, rescale=True, format='norm', firstdims=None, **kwargs):

        # Unstack: level 0
        stacks = npy.split(sdata, self.splits, axis=0)

        # Ustack: level 1
        unstacks= [self[i].unstack(pdata, rescale=rescale, format=format,
            firstdims=firstdims, **kwargs) for i, pdata in enumerate(stacks)]

        return self.unmap(unstacks)

    def project_model(self, var, checkid=True):
        """Project model variables to observations positions


        .. warning:: This methods projects your variable on ALL platforms,
            even if the variable has not observation equivalent.
            It will emit a warning if ``checkid`` is True and the id of
            the variable is not known from the platform.
            See :meth:`~sonat.ens.Ensemble.project_on_obs` for a valid
            projection approach.

        """
        return self.unmap([obs.project_model(var, checkid=checkid)
                           for obs in self])

    def assert_compatible_with_ens(self, ens, syncnorms=True):
        """Assert that an :class:`~sonat.obs.Ensemble` current instance is compatible
        with the current :class:`ObsManager` instance

        It checks that observed variable are provided by the ensemble.
        It optinonally synchonise norms between model and observations.
        """
        ens.assert_compatible_with_obs(self, syncnorms=syncnorms)

    def plot(self, variables=None, input_mode='names',
             full3d=True, full2d=True,
             lon_bounds_margin=.1, lat_bounds_margin=.1,
             zonal_sections=None, merid_sections=None, horiz_sections=None,
             level=None, lon=None, lat=None,
             sync_vminmax=True,
             color=None, marker=None, legend=True,
             color_cycle='bgrcmy', marker_cycle='o^s<>*',
             reset_cache=True, fig=None, savefig=True, close=True,
             figpat='sonat.obs.{var_name}_{slice_type}_{slice_loc}.png',
             subst={}, **kwargs):
        """Plot the locations of all platforms

        See :meth:`NcObsPlatform.plot` for more arguments.


        Parameters
        ----------
        lon: tuple of floats
            Longitude bounds
        lat: tuple of floats
            Latitude bounds
        level: tuple of floats
            Level bounds
        lon_bounds_margin: float
            Fraction of longitude bounds to add as margin when bounds
            are not explicitely specified
        lat_bounds_margin: float
            Fraction of latitude bounds to add as margin when bounds
            are not explicitely specified
        level_bounds_margin: float
            Fraction of level bounds to add as margin when bounds
            are not explicitely specified
        sync_vminmax: bool
            Sync min and max between different plots of the same variables
        """
        self.verbose('Plotting observations locations')

        # Input mode
        valid_input_modes = ['names', 'arrays']
        if input_mode not in valid_input_modes:
            raise SONATError('input_mode must be one of: ' +
                ', '.join(valid_input_modes))
        if input_mode=='arrays':
            variables = self.remap(variables)

        # Cache
        if reset_cache:
            del self.plot_cache

        # Cycler
        if color is None and marker is None:
            if len(color_cycle)!=len(marker_cycle) or len(color_cycle)<len(self):
                cyc = (cycler(color=color_cycle) * cycler(marker=marker_cycle))()
            else:
                cyc = cycler(color=color_cycle, marker=marker_cycle)()
        elif color is None:
            cyc = cycler(color=color_cycle)()
        else:
            cyc = cycler(marker=marker_cycle)()

        # Bounds
        if lon is None and (full2d or full3d or zonal_sections or horiz_sections):
            lon = self.get_lon(margin=lon_bounds_margin)
        if lat is None and (full2d or full3d or merid_sections or horiz_sections):
            lat = self.get_lat(margin=lat_bounds_margin)
        if (level is None and (not self.has_bottom or bathy is not None) and
            (full3d or zonal_sections or merid_sections)):
            level = (self.get_level()[0], 0)

        # Loop on platforms
        figs = {}
        for ip, obs in enumerate(self):

            # Check var names
            if input_mode=='names':
                myvariables = obs.get_valid_varnames(variables)
                if not myvariables:
                    if variables:
                        continue
                    else: # nothing requested, no locations
                        myvariables = ['locations']
            else:
                myvariables = variables[ip]

#            # Loops status
#            isfirst = ip==0
#            islast = ip==len(self)-1

            # Sync caching
            obs.plot_cache = self.plot_cache

            # Symbols
            kwargs.update(color = dicttree_get(color, obs.name),
                marker = dicttree_get(marker, obs.name))
            if kwargs['color'] is None or kwargs['marker'] is None:
                c = cyc.next()
                if kwargs['color'] is None:
                    kwargs['color'] = c['color']
                if kwargs['marker'] is None:
                    kwargs['marker'] = c['marker']

            # Plot
            obs.plot(variables=myvariables, reset_cache=False,
                     full2d=full2d, full3d=full3d,
                     lon=lon, lat=lat, level=level,
                     fig=fig, figpat=figpat,
                     savefig=False, close=False,
                     legend='cache' if legend else False,
                     sync_vminmax=sync_vminmax,
                     **kwargs)

        # Plot all pending legends
        if legend:
            self.add_cached_plots_legend()

        # Save all plots
        if savefig:
            kw = locals().copy()
            del kw['self']
#            kw['subst'] = subst
            kw.update(subst)
            return self.save_cached_plots(**kw)
        return {}



def load_obs(ncfiles, varnames=None, lon=None, lat=None,
        logger=None):
    """Quickly load all observations with one file per platform"""
    # Logger
    if logger is None:
        logger = get_logger()

    # Load platforms
    if isinstance(ncfiles, basestring):
        ncfiles = [ncfiles]
    obsplats = [NcObsPlatform(ncfile, varnames=varnames, lon=lon, lat=lat,
        logger=logger) for ncfile in ncfiles]

    # Init ObsManager
    return ObsManager(obsplats, logger=logger)


register_obs_platform(NcObsPlatform)
