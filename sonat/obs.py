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
import cdms2
import MV2
import numpy as N
from vcmq import (grid2xy, regrid2d, ncget_lon, ncget_lat,
    ncget_time, ncget_level, ArgList, ncfind_obj, itv_intersect, intersect,
    MV2_axisConcatenate, transect, create_axis)

from .__init__ import sonat_warn, SONATError, BOTTOM_VARNAMES
from .misc import xycompress, _Base_, _XYT_, check_variables
from .stack import Stacker, _StackerMapIO_

npy = N

RE_PERTDIR_MATCH = re.compile(r'[+\-][xy]$').match

class NcObsPlatform(Stacker, _XYT_):
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

    def __init__(self, ncfile, varnames=None, logger=None, norms=None,
            time=None, lon=None, lat=None, levels=None,
            singlevar=False, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Init parameters
        self.ncfile = ncfile
        self.varnames = ArgList(varnames).get() if varnames else varnames
        self.time = time
        self.lon = lon
        self.lat = lat
        self.levels = levels
        self.errors = OrderedDict()
        self._orig = {}

        # Load positions and variables
        self.load(singlevar=singlevar)

        # Init stacker
        Stacker.__init__(self, self.errors.values(), logger=False, means=False,
            norms=norms)


    def load(self, singlevar=False, **kwargs):
        """Read mobility, errors, depth and time in open netcdf file"""
        # Open
        f = cdms2.open(self.ncfile)

        # List of error variables
        if self.varnames is None:
            self.varnames = [vname.rstrip(self.nc_error_suffix)
                for vname in f.listvariables()
                    if vname.endswith(self.nc_error_suffix)]
        else:
            self.varnames = [vname.rstrip(self.nc_error_suffix) for vname in self.varnames]
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
        if grid and len(grid.shape)==2: # Structured grid

            self.pshape = 'gridded'
            kwread = dict(time=self.time, lat=self.lat, lon=self.lon)
            if '-' in order:
                raise SONATError("There are unkown dimensions in your gridded variable. "
                    "Current order: "+order)
            self.axes1d = order[:-2]


        else: # Unstructured grid

            mask = N.ma.nomask #zeros(fs.shape[-1], '?')
            self.pshape = 'xy'
            self.axes1d = []
            kwread = {}

            # Lon/lat selection
            if grid:
                kwread.update(lat=self.lat, lon=self.lon)
            else:
                lons = ncget_lon(f)
                if lons is None:
                    raise SONATError("Longitudes not found")
                lats = ncget_lat(f)
                if lats is None:
                    raise SONATError("Latitudes not found")
                lons = lons.asma()
                lats = lats.asma()
                if self.lon:
                    mask |= lons < self.lon[0]
                    mask |= lons > self.lon[1]
                if self.lat:
                    mask |= lats < self.lat[0]
                    mask |= lats > self.lat[1]
            mask = mask.filled(False)
            order = order[:-1]

            # Vertical dimension
            if 'z' in order:
                order.remove('z')
                self.depths = fs.getLevel().clone()
                self.axes1d.append(self.depths)
            else:
                self.depths = ncget_level(f)
                if self.depths is not None:
                    raise SONATError("Scattered Z dimension not yet supported")
                    self.pshape = self.pshape+'z'
                elif hasattr(f, 'depth'): # depth from file attribute
                    self.depths = f.depth
                    if not isinstance(self.depths, basestring):
                        raise SONATError("Depth must be a string if specified as an attribute")
                elif fs.id in BOTTOM_VARNAMES: # 2D: bottom
                    self.depths = 'bottom'
                else: # 2D: surf
                    self.depths = 'surf'

            # Time selection
            if 't' in order: # 1D axis
                kwread.update(time=time)
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
        for vname in self.varnames:
            self.errors[vname] = f(vname + self.nc_error_suffix, **kwread)
        self._sample = sample = self.errors[vname]
        # - lon/lat
        if grid:
            self.lons = sample.getLongitude()[:]
            self.lats = sample.getLatitude()[:]
        else:
            self.lons = lons
            self.lats = lats
        # - mobility
        if self.nc_mobility_name in f.listvariables():
            self.mobility = f(self.nc_mobility_name, **kwread)
        elif self.pshape=='gridded': # don't move grid for the moment
            self.mobility = 0.
        else:
            sonat_warn('No mobility variable found in obs file. '
                'Setting it to 1 everywhere.')
            shape = grid.shape if grid else lons.shape
            self.mobility = MV2.ones(shape, 'i')
        self.mobile = (N.ma.asarray(self.mobility)==1).all()

        # Compression to remove masked point
        if self.pshape != 'gridded' and mask.any():
            valid = ~mask
            self.lons = xycompress(valid, self.lons)
            self.lats = xycompress(valid, self.lats)
            for vname, var in self.errors.items():
                self.errors[vname] = xycompress(valid, var)
            if self.times is not None and self.times not in self.axes1d:
                self.times = xycompress(valid, self.times)
            if (self.depths is not None and self.depths not in self.axes1d and
                    not isinstance(self.depths, basestring)):
                self.depths = xycompress(valid, self.depths)

        f.close()



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

    def get_error(vname):
        return self.errors[vname]

    def set_named_norms(self, *anorms, **knorms):
        """Set norms by variable names

        Note
        ----
        The :class:`~sonat.pack.Packer` instance of  allvariables that were
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


    def interp_model(self, var):
        """Interpolate model variables to observations positions"""
        # List of variables
        al = ArgList(var)

        # Loop on variables
        out = []
        for i, var in enumerate(al.get()):

            # Check id
            if var.id not in self.varnames:
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
            for axis in self.axes1d:
                var = regrid1d(var, axis, method='linear')

            out.append(var)

        return al.put(out)

#    def plot(self):
#        pass

class ObsManager(_Base_, _StackerMapIO_, _XYT_):
    """Class to manage several observation platform instances"""

    def __init__(self, input, logger=None, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Load platforms
        obsplats = _StackerMapIO_._load_input_(self, input)
        self.obsplats = []
        for obsplat in obsplats:
            if not isinstance(obsplat, NcObsPlatform):
                raise SONATError('ObsManager must be initialised with a single or '
                    'list of NcObsPlatform instances')
            self.obsplats.append(obsplat)

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

    @property
    def varnames(self):
        vv = []
        for obs in self:
            vv.extend(obs.varnames)
        return list(set(vv))

    @property
    def lons(self):
        if not hasattr(self, '_lons'):
            self._lons = N.concatenate([obs.lons[:].ravel() for obs in self
                if obs.lons is not None])
        return self._lons

    @property
    def lats(self):
        if not hasattr(self, '_lats'):
            self._lats = N.concatenate([obs.lats[:].ravel() for obs in self
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
            if varname in self.varnames:
                if not self._depths[varname]:
                    self._depths[varname] = None
                else:
                    self._depths[varname] = tuple(self._depths[varname])

        return self._depths

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

    def interp_model(self, var):
        """Interpolate model variables to observations positions"""
        return self.unmap([obs.interp_model(var) for obs in self])


