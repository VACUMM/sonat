"""
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
import re
from vcmq import (cdms2, MV2, grid2xy, regrid2d, ncget_lon, ncget_lat,
    ncget_lat, ncget_dep, ArgList, ncfind_obj, itv_intersect, intersect)

from .__init__ import _Base_
from stack import Stacker


RE_PERTDIR_MATCH = re.compile(r'[+\-][xy]$').match

class NcObsPlatform(Stacker):
    """Generic observation platform class


    Attributes
    ----------
    pshape: string
        Plateform shape
    """

    nc_error_suffix = '_error'
    nc_flag_name = 'flag'

    def __init__(self, ncfile, ncvars=None, logger=None, norms=None,
            time=None, lon=None, lat=None, levels=None, pert=1e-2, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Init parameters
        self.ncfile = ncfile
        self.ncvars = ncvars
        self.time = time
        self.lon = lon
        self.lat = lat
        self.levels = levels
        self.errors = OrderedDict()
        self.pert = pert
        self._orig = {}

        # Load positions and variables
        self.load()

        # Init stacker
        Stacker.__init__(self.errors.values(), logger=False, mean=False,
            norms=norms)


    def load(self, **kwargs):
        """Read flag, errors, depth and time in opened netcdf file"""
        # Open
        f = cdms2.open(self.ncfile)

        # List of error variables
        if self.ncvars is None:
            self.ncvars = [vname for vname in f.listvariables()
                if vname.endswith(nc_error_suffix)]
        else:
            self.ncvars = [vname.lstrip(nc_error_suffix) for vname in self.ncvars]
        if not self.ncvars:
            raise PyARMError(('No valid error variable with suffix "{self.nc_error_suffix}"'
                ' in file: {self.ncfile}').format(self))

        # Reference variable
        fs = f[self.ncvars[0] + nc_error_suffix]

        # Inspect
        grid = fs.getGrid()
        order = fs.getOrder()
        if grid and len(grid.shape)==2: # Structured grid

            self.pshape = 'gridded'
            kwread = dict(time=self.time, lat=self.lat, lon=self.lon)
            if '-' in self.order:
                raise PyARMError("There are unkown dimensions in your gridded variable. "
                    "Current order: "+self.order)
            self.axes1d = self.order[:-2]

        else: # Unstructured grid

            mask = N.ma.nomask #zeros(fs.shape[-1], '?')
            self.pshape = 'xy'
            self.axes1d = []

            # Lon/lat selection
            if grid:
                kwread = dict(lat=self.lat, lon=self.lon)
#                lons = grid.getLongitude()
#                lats = grid.getLatitude()
            else:
                lons = ncget_lon(f)
                if lons is None:
                    raise PyARMError("Longitudes not found")
                lons = ncread_lon(f, id=lons)
                lats = ncget_lat(f)
                if lats is None:
                    raise PyARMError("Latitudes not found")
                if self.lon:
                    mask |= lons < self.lons[0]
                if self.lat:
                    mask |= lats < self.lats[1]
            order = order[:-1]

            # Vertical dimension
            if 'z' in self.order:
                order.remove('z')
                self.depths = fs.getLevel().clone()
                self.axes1d.append(self.depths)
            else:
                self.depths = ncget_dep(f)
                if self.depths is not None:
                    raise PyARMError("Scattered Z dimension not yet supported")
                    self.pshape = self.pshape+'z'
                elif hasattr(fs, 'depth'): # depth from attribute
                    self.depths = fs.depth
                    if not isinstance(self.depths, basestring):
                        raise PyARMError("Depth must be a string if specified as an attribute")

            # Time selection
            if 't' in self.order: # 1D axis
                kwread = dict(time=time)
                order.remove('t')
                self.times = fs.getTime().clone()
                self.axes1d.append(self.times)
            else: # Aux axis
                self.times = ncget_time(f)
                if self.times is not None:
                    if self.time:
                        times = create_time(times[:], times.units)[:]
                        mask |= times < reltime(self.time[0], times.units)
                        mask |= times > reltime(self.time[1], times.units)
                self.pshape = self.pshape + 't'

            # Check mask
            if mask.all():
                PyARMError("All your observation data are masked")

            # Check remaining dims
            if order:
                PyARMError("There are unkown dimensions in your scattered obs variable")


        # Read
        # - errors
        for vname in self.ncvars:
            self.errors[vname] = f(vname + nc_error_suffix, **kwread)
        sample = self.errors[vname]
        # - lon/lat
        if grid:
            self.lons = sample.getLongitude()[:]
            self.lats = sample.getLatitude()[:]
        else:
            self.lons = lons
            self.lats = lats
        # - flag
        if self.nc_flag_name in f.listvariables():
            self.flag = f(self.nc_flag_name, **kwread)
        else:
            pyarm_warn('No flag variable found in obs file. '
                'Setting it to 1 everywhere.')
            shape = grid.shape if grid else lons.shape
            self.flag = MV2.ones(shape, 'i')

        # Compression to remove masked point
        if self.pshape != 'gridded' and mask.any():
            valid = ~mask
            self.lons = xycompress(valid, self.lons)
            self.lats = xycompress(valid, self.lats)
            for vname, var in self.errors.item():
                self.errors[vname] = xycompress(valid, var)
            if self.times is not None and self.times not in self.axes1d:
                self.times = xycompress(valid, self.times)
            if self.depths is not None and self.depths not in self.axes1d:
                self.depths = xycompress(valid, self.depths)



    @property
    def ndim(self):
        return self.flag.ndim

    @property
    def shape(self):
        return self.flag.shape

    @property
    def grid(self):
        return self.flag.getGrid()

    @property
    def ctimes(self):
        if not hasattr(self, '_ctimes'):
            if self.times is None:
                self._ctimes = None
            else:
                self._ctimes = comptime(self.times)
        return self._ctimes

    def get_seldict(self, axes='xyt', xybounds='cce', tbounds='cce'):
        sel = {}
        if 'x' in axes:
            sel['lon'] = (self.lons.min(), self.lons.max(), xybounds)
        if 'y' in axes:
            sel['lat'] = (self.lats.min(), self.lats.max(), xybounds)
        if 't' in axes and self.ctimes:
            sel['time'] = (self.ctimes.min(), self.ctimes.max(), tbounds)
        return sel

    def intersects(self, lon=None, lat=None, time=None):
        """Does the observations intersects specified intervals"""
        sdict = self.get_seldict()
        if lon is not None and not intersect(lon, sdict['lon']):
            return False
        if lat is not None and not intersect(lat, sdict['lat']):
            return False
        if (time is not None and self.times is not None and
                not intersect(time, sdict['time'])):
            return False
        return True


    def get_error(vname):
        return self.errors[vname]

    def activate_xy_pert(self, direction, index, mres):
        """

        Parameters
        ----------
        direction: string
            Direction of perturbation in the form {+|-}{x|y}"
        mgrid: cdms2 grid
            Model grid
        """
        if self.rshape=='gridded':
            return

        if direction:

            # Parse
            direction = str(direction).lower()
            if not RE_PERTDIR_MATCH(direction):
                raise PyARMError("The direction of perturbation argument must"
                    " be of the form: '{+|-}{x|y}'")
            isx = direction[1]=='x'
            sign = 1 if direction[0]=='-' else -1

            # X and Y perturbation values
            if not N.isscalar(mres): # grid
                mres = min(resol(mres, meters=False))
            base_pert = self.pert * mres * sign
            if isx:
                pert = base_pert / N.radians(self.lats).mean()
                cname = 'lons'
            else:
                pert = base_pert
                cname = 'lats'
            coords = getattr(self, cname) # array of coordinates

            # Reset
            self.deactivate_pert()

            # Save original coordinates values
            self._orig[cname] = coords.copy()

            # Change coordinates in place
            if coord.ndim==2:
                index = N.unreavel_index(coords.shape)
            coords[index] += pert

        elif self._orig: # back to orig

            self.deactivate_pert()


    def deactivate_xy_pert(self):

        for cname in 'lons', 'lats':
            if cname in self._orig:
                getattr(sel, cname)[:] = self._orig[cname]
                del self._orig[cname]

    @property
    def xy_pert_indices(self):
        if hasattr(self, '_xy_pert_indices'):
            return self.xy_pert_indices
        self.xy_pert_indices = N.where(self.flag>=1)[0]
        return self.xy_pert_indices


    def interp_model(self, var):

        # List of variables
        al = ArgList(var)

        # Loop on variables
        out = []
        for var in al.get():

            # Order
            order = var.getOrder()

            # Gridded and scattered
            if self.rtype=='gridded':
                var = regrid2d(var, self.grid, method='bilinear')
            else:
                kw = {}
                if 't' in self.pshape:
                    if 't' not in order:
                        raise PyARMError("Model variables must have a time axis")
                    kw['times'] = self.times
                if 'z' in self.pshape:
                    if 'z' not in order:
                        raise PyARMError("Model variables must have a depth axis")
                    kw['depths'] = self.depths
                var = transect(var, lons=lons, lats=lats, **kw)

            # Auxilary axes
            for axis in self.axes1d:
                var = regrid1d(var, axis, method='linear')

            out.append(var)

        return al.put(out)

    def plot(self):
        pass

def xycompress(valid, vari):
    """Keep valid spatial points"""
    # Init
    nv = valid.sum()
    ax = vari.getAxis(-1)
    vari = vari[:nv].clone()

    # Fill
    vari.getAxis(-1)[:] = N.compress(valid, ax[:])
    vari[:] = N.compress(valid, vari.asma(), axis=-1)

    return vari

class ObsManager(_Base_):

    def __init__(self, obs, logger=None, **kawargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        self.stacked_data = npy.asfortranarray(
            npy.concatenate([p.packed_data for p in self.packers], axis=0))
        self.splits = npy.cumsum([p.packed_data.shape[0] for p in self.packers[:-1]])
