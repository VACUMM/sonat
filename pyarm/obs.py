"""
"""

from vcmq import cdms2, MV2, grid2xy, regrid2d

from .__init__ import _Base_

class ObsPlatform(_Base_):

    nc_error_prefix = 'error_'
    nc_flag = 'flag'

    def __init__(self, ncfile, ncvars=None, logger=None,
            lon=None, lat=None, levels=None, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Init parameters
        self.ncfile = ncfile
        self.ncvars = ncvars
        self.lon = lon
        self.lat = lat
        self.levels = levels
        self.errors = {}

        # Load positions and variables
        self.load()

    def load(self):
        pass

    def _load_variables_(self, f, **kwargs):
        """Read flag and errors in opened netcdf file"""

        self.flag = f(self.nc_flag, **kwargs)

        if self.ncvars is None:
            self.ncvars = [vname for vname in f.listvariables()
                if vname.startswith(nc_error_prefix)]
        self.ncvars = [vname.lstrip(nc_error_prefix) for vname in self.ncvars]
        if not self.ncvars:
            self.warning('No valid error variable in: '+self.ncfile)
        else:
            for vname in self.ncvars:
                self.errors[vname] = f(nc_error_prefix+vname, **kwargs)
        
    def get_seldict(self, axes='xy', bounds='cce'):
        sel = {}
        if 'x' in axes:
            sel['lon'] = (self.lons.min(), self.lons.max(), xyb)
        if 'y' in axes:
            sel['lat'] = (self.lats.min(), self.lats.max(), xyb)
        return sel


    def get_error(vname):
        return self.errors[vname]

    def plot(self):
        pass

    def stack(self, var):
        if isinstance(var, basestring):
            var = self[var]
        pass

    def unstack(self, arr):
        pass

    def get_xy_perturbation(self):
        pass

    def interp_model(self, var):
        pass

class Unstruct(object):
    
    nc_lon = 'lon'
    nc_lat = 'lat'

    def load(self):

        logger.debug('Loading file: ' + self.ncfile)

        # Read
        f = cdms2.open(self.ncfile)
        self.lons = f(self.nc_lon, raw=True)
        self.lats = f(self.nc_lat, raw=True)
        self._load_variables_(f)
        f.close()

        # Mask
        if self.lon is not None or self.lat is not None:
            mask = self.lons.mask | self.lats.mask
            if self.lon is not None:
                mask |= self.lons < self.lon[0]
            if self.lat is not None:
                mask |= self.lats < self.lat[0]
            for var in self.errors.values() + [self.flag]:
                var[:] = MV2.masked_where(mask, var, copy=False)

        

class SurfacePoints(ObsPlatform, Unstruct):

    

        