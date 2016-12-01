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
import re
from vcmq import (cdms2, MV2, grid2xy, regrid2d, ncget_lon, ncget_lat,
    ncget_lat, ncget_dep, ArgList, ncfind_obj, itv_intersect, intersect,
    MV2_axisConcatenate)

from .misc import xycompress, _Base_, _XYT_
from .stack import Stacker


RE_PERTDIR_MATCH = re.compile(r'[+\-][xy]$').match

class NcObsPlatform(Stacker, _XYT_):
    """Generic observation platform class


    Attributes
    ----------
    pshape: string
        Plateform shape
    """

    nc_error_suffix = '_error'
    nc_flag_name = 'flag'

    def __init__(self, ncfile, ncvars=None, logger=None, norms=None,
            time=None, lon=None, lat=None, levels=None, pert=1e-2,
            singlevar=False, **kwargs):

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
        self.load(singlevar=singlevar)

        # Init stacker
        Stacker.__init__(self.errors.values(), logger=False, means=False,
            norms=norms)


    def load(self, singlevar=False, **kwargs):
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
            raise SONATError(('No valid error variable with suffix "{self.nc_error_suffix}"'
                ' in file: {self.ncfile}').format(self))
        if singlevar:
            self.ncvars = self.ncvars[:1]

        # Reference variable
        fs = f[self.ncvars[0] + nc_error_suffix]

        # Inspect
        grid = fs.getGrid()
        order = fs.getOrder()
        if grid and len(grid.shape)==2: # Structured grid

            self.pshape = 'gridded'
            kwread = dict(time=self.time, lat=self.lat, lon=self.lon)
            if '-' in self.order:
                raise SONATError("There are unkown dimensions in your gridded variable. "
                    "Current order: "+self.order)
            self.axes1d = self.order[:-2]

        else: # Unstructured grid

            mask = N.ma.nomask #zeros(fs.shape[-1], '?')
            self.pshape = 'xy'
            self.axes1d = []

            # Lon/lat selection
            if grid:
                kwread = dict(lat=self.lat, lon=self.lon)
            else:
                lons = ncget_lon(f)
                if lons is None:
                    raise SONATError("Longitudes not found")
                lons = ncread_lon(f, id=lons)
                lats = ncget_lat(f)
                if lats is None:
                    raise SONATError("Latitudes not found")
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
                    raise SONATError("Scattered Z dimension not yet supported")
                    self.pshape = self.pshape+'z'
                elif hasattr(f, 'depth'): # depth from file attribute
                    self.depths = f.depth
                    if not isinstance(self.depths, basestring):
                        raise SONATError("Depth must be a string if specified as an attribute")

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
                SONATError("All your observation data are masked")

            # Check remaining dims
            if order:
                SONATError("There are unkown dimensions in your scattered obs variable")


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
            sonat_warn('No flag variable found in obs file. '
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
                raise SONATError("The direction of perturbation argument must"
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
        """Interpolate model variables to observations positions"""
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
                        raise SONATError("Model variables must have a time axis")
                    kw['times'] = self.times
                if 'z' in self.pshape:
                    if 'z' not in order:
                        raise SONATError("Model variables must have a depth axis")
                    kw['depths'] = self.depths
                var = transect(var, lons=lons, lats=lats, **kw)

            # Auxilary axes
            for axis in self.axes1d:
                var = regrid1d(var, axis, method='linear')

            out.append(var)

        return al.put(out)

    def plot(self):
        pass

class ObsManager(_Base_):

    def __init__(self, input, logger=None, **kawargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Load platforms
        obsplats = _MapIO_._load_input_(input)
        self.obsplats = []
        for obsplat in obsplats:
            if not isinstance(obsplat, NcObsPlatform):
                raise SONATError('ObsManager must be initialised with a single or '
                    'list of NcObsPlatform instances')
                self.obsplats.append(obsplat)

        # Stack stacked data
        self.stacked_data = npy.asfortranarray(
            npy.concatenate([o.stacked_data for o in self.obsplats], axis=0))
        self.splits = npy.cumsum([o.stcked_data.shape[0] for o in self.obsplats[:-1]])


    def __len__(self):
        return len(self._obsplats)

    def __getitem__(self, key):
        return self.obsplats[key]

    @property
    def ncvars(self):
        vv = []
        for obs in self:
            vv.extend(obs.ncvars)
        return list(set(vv))

    @property
    def lons(self):
        if not hasattr(self, '_lons'):
            self._lons = N.concatenate(obs.lons[:].ravel() for obs in self
                if obs.lons is not None)
        return self._lons

    @property
    def lats(self):
        if not hasattr(self, '_lats'):
            self._lats = N.concatenate(obs.lats[:].ravel() for obs in self
                if obs.lats is not None)
        return self._lats

    @property
    def times(self):
        if not hasattr(self, '_lats'):
            self._lats = MV2_axisConcatenate(obs.times for obs in self
                if obs.times is not None)
        return self._times

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

    def interp_model(self, var):
        """Interpolate model variables to observations positions"""
        return sel.unmap([obs.interp_model(var) for obs in self])

