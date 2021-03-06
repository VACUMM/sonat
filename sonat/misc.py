"""Misc utilities"""
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
import string
from string import Formatter
from glob import has_magic, glob
from collections import OrderedDict
import numpy as N
import cdms2
import MV2
from genutil import minmax as minmax
from vcmq import (ncget_time, pat2freq, lindates, adatetime,
    comptime, add_time, pat2glob, are_same_units, indices2slices,
    kwfilter, GENERIC_VAR_NAMES, DS, set_atts, format_var,
    match_known_var, ArgList, create_lon, regrid1d, grid2xy, create_lat,
    create_time, create_dep, create_axis, cp_atts, isaxis,
    dicttree_set, dicttree_get, intersect)

from .__init__ import sonat_warn, SONATError, get_logger

def scan_format_string(format_string):
    """Scan a format string using :class:`string.Formatter`

    Parameters
    ----------
    format_string: string

    Return
    ------
    dict:
        Positional and keyword fields with their properties in a dict
    dict:
        Dictionary of properties with the follwing keys:

        - positional: list of positional field keys
        - keyword: list of keyword field keys
        - with_time: list of keys that have a date pattern format
    """
    fields = {}
    props = {'with_time': [], 'positional':[], 'keyword':[]}
    f = Formatter()
    for literal_text, field_name, format_spec, conversion in f.parse(format_string):
        if field_name is None:
            continue
        first, _ = type(field_name)._formatter_field_name_split(field_name)

        scan = dict(field_name=field_name, format_spec=format_spec, conversion=conversion)
        fields[first] = scan

        if isinstance(first, (int, long)):
            props['positional'].append(first)
        else:
            props['keyword'].append(first)
        if '%' in format_spec:
            props['with_time'].append(first)

    return fields, props


class DatePat2GlobFormatter(string.Formatter):
    def _vformat(self, format_string, args, kwargs, used_args, recursion_depth):
        if recursion_depth < 0:
            raise ValueError('Max string recursion exceeded')
        result = []
        for literal_text, field_name, format_spec, conversion in \
                self.parse(format_string):

            # output the literal text
            if literal_text:
                result.append(literal_text)

            # if there's a field, output it
            if field_name is not None:
                # this is some markup, find the object and do
                #  the formatting
#                print 'list', literal_text, '|', field_name, '|', format_spec, '|', conversion
                # given the field_name, find the object it references
                #  and the argument it came from
                try:
                    obj, arg_used = self.get_field(field_name, args, kwargs)
                except:
                    if format_spec is not None and '%' in format_spec:
                        result.append(pat2glob(format_spec))
                        continue
                used_args.add(arg_used)

                # do any conversion on the resulting object
                obj = self.convert_field(obj, conversion)

                # expand the format spec, if needed
                format_spec = self._vformat(format_spec, args, kwargs,
                                            used_args, recursion_depth-1)

                # format the object and append to the result
                result.append(self.format_field(obj, format_spec))

        return ''.join(result)



def list_files_from_pattern(ncpat, time=None, dtfile=None, sort=True, **subst):
    """List files possibly with glob and date patterns

    Parameters
    ----------
    ncpat: string
        File name with date patterns
    time: tuple, None
        Date interval
    dtfile: tuple, None
        Time step between two files like ``(10, 'days')``.
        This time step is assumed to be constant across files.
    sort: bool
        Sort after listing?
    \**subst: dict
        Use for substitution in ``ncpat``.
    """

    # List all files
    if isinstance(ncpat, list): # A list of file

        files = []
        for filepat in ncpat:
            files.extend(list_files_from_pattern(filepat, time=time, dtfile=dtfile, **subst))

    else: # A single string

        with_magic = has_magic(ncpat)

        scan_fields, scan_props = scan_format_string(ncpat)
        if scan_props['with_time']: # With time pattern

            # With time
            if time is None: # without

                sonat_warn("You should better provide a time interval "
                    "with a date pattern in file name")
                ncfile = DatePat2GlobFormatter().format(ncpat, **subst)
                files = glob(ncfile)

            else: # with

                # Guess pattern and frequency
                date_format = scan_fields[scan_props['with_time'][0]]['format_spec']
                freq = pat2freq(date_format)
                if dtfile is None:
                    dtfile = 1, freq
                    sonat_warn('Time steps between files not explicitly specified. '
                        'Set to {}. You may miss first files!'.format(dtfile))
                elif not isinstance(dtfile, tuple):
                    dtfile = dtfile, freq

                # Generate dates or glob patterns
                files = []
                ct0 = add_time(time[0], -dtfile[0], dtfile[1])
                ct1 = time[-1]
                for date in lindates(ct0, ct1, 1, dtfile[1]):
                    date = adatetime(date)
                    ss = subst.copy()
                    ss['date'] = date
                    ncfile = ncpat.format(**ss)
                    if with_magic:
                        files.extend(glob(ncfile))
                    elif os.path.exists(ncfile):
                        files.append(ncfile)

        elif has_magic(ncpat): # Just glob pattern

                files = glob(ncpat)

        else: # Just a file

                files = [ncpat]

    # Check existence
    files = filter(os.path.exists, files)

    # Unique
    files = list(set(files))

    # Sort
    if sort:
        files.sort(key=sort if callable(sort) else None)

    return files

def ncfiles_time_indices(ncfiles, dates, getinfo=False, asslices=False):
    """Get time indices corresponding to each dates for each files

    .. warning:: The time axis is read only for the first two files.
        All files must obey the following rules:

        - Be in chronological order.
        - Have the same time step.
        - Have the same number of time steps.

        Dates are sorted chronologically before processing.

    """

    # Dates
    dates = comptime(dates)
    dates.sort()

    # Select needed files
    if not ncfiles:
        return []
    ncfdict = OrderedDict()
    duplicate_dates = []
    for i, ncfile in enumerate(ncfiles):

        # Get file times
        if i<2: # Read time

            taxis = ncget_time(ncfile, ro=True)
            if taxis is None:
                SONATError("Can't read time axis in file: " + ncfile)
#            ctimes = taxis.asComponentTime()

            if i==0: # Reference info

                t0 = taxis[0]
                taxis0 = taxis.clone()
                tunits = taxis0.units

            else: # Get file time interval

                if not are_same_units(taxis.units, tunits):
                    taxis.toRelativeTime(tunits)
                dt = taxis[0] - t0

        else: # Generate time

            taxis = taxis0.clone()
            taxis[:] += dt * i

        # Loop on dates
        for date in list(dates): #enumerate(list(dates)):

            # Get index
            ijk = taxis.mapIntervalExt((date, date), 'cob')

            # Checks
            if ijk is None: # Date not in file
                continue
            it = ijk[0]
            if ncfile not in ncfdict: # Init
                ncfdict[ncfile] = [it]
            elif it in ncfdict[ncfile]: # Time step already used
                duplicate_dates.append(date)
            else:
                ncfdict[ncfile].append(it)
            dates.remove(date)

        # As sclices
        if asslices and ncfile in ncfdict:
            ncfdict[ncfile] = indices2slices(ncfdict[ncfile])

    if not getinfo:
        return ncfdict
    return ncfdict, {'missed':dates, 'duplicates':duplicate_dates}

def asma(*args):
    """Return pure numpy or numpy.ma arrays"""
    if not args: return
#    single = len(args)==1
    out = []
    for arg in args:
        if cdms2.isVariable(arg):
            arg = arg.asma()
        out.append(arg)
    if len(args)==1:
        return out[0]
    return out

class NcReader(object):
    """A generic interface to open and read a netcdf file

    Examples
    --------

    >>> f = NcReader(ncfile, 'mars3d') # VACUMM Mars3D reader
    >>> u = f('u', lon=(10, 12))    # generic var name
    >>> u = f('+sst', lon=(10, 12)) # netcdf var name
    >>> f.close()

    >>> f = NcReader(ncfile, 'cdms2') # direct cdms2 reader
    >>> u = f('u', lon=(10, 12)) # netcdf var name
    >>> f.close()

    >>> f = NcReader(ncfile, 'netcdf') # python netcdf4 reader
    >>> u = f('u', slice(0, 10)) # <=> f['u'][:10]
    >>> f.close()

    """

    def __init__(self, ncfile, readertype='generic',  **kwargs):
        self.ncfile = ncfile
        if not isinstance(readertype, basestring):
            self.f = readertype(ncfile,  **kwargs)
            self.type = 'callable'
        elif readertype=='cdms2':
            self.f = cdms2.open(ncfile,  **kwargs)
            self.type = readertype
        elif readertype=='netcdf4':
            import netcdf4
            self.f = netcdf4.Dataset(ncfile,  **kwargs)
            self.type = readertype
        else:
            self.f = DS(ncfile, readertype, **kwargs)
            self.type = 'dataset'


    def __call__(self, vname, *args, **kwargs):

        if self.type=='netcdf4':

            args = tuple(filter(lambda x: isinstance(x, slice), args))

            return self.f[vname][args]

        else:

            return self.f(vname, *args, **kwargs)

    def __getitem__(self, vname):
        return self.f[vname]

    def get_variables(self):
        if self.type=='netcdf4' or self.type=='cdms2':
            return self.variables.keys()
        if self.type=='dataset':
            return self.f.dataset[0].variables.keys()
        raise SONATError('Method yet not implemented for read type "{}"'.format(self.type))

    def close(self):
        self.f.close()


def xycompress(valid, vari, **atts):
    """Keep valid spatial points

    Parameters
    ----------
    valid: 1D or 2D bool array
        2D array for data on structured grid
    vari: array
        The variable to compress
    """

    # Slice for extra dimensions
    pslice = (slice(None), ) * (vari.ndim - valid.ndim)

    if cdms2.isVariable(vari):

        nv = valid.sum()


        if valid.ndim==2:


            # Init
            assert valid.ndim == 2, 'Valid must be a 2D array'
            varo = vari[pslice + (0, slice(0, nv))].clone()
            ax = create_axis((nv, ), id='point', long_name='Spatial points')
            varo.setAxis(-1, ax)
            varo.setGrid(None)

            # Fill
            varo[:] = vari.asma()[pslice + (valid, )]

        else:

            # Init
            ax = vari.getAxis(-1)
            varo = vari[pslice + (slice(0, nv), )].clone()

            # Fill
            varo.getAxis(-1)[:] = N.compress(valid, ax[:])
            varo[:] = N.ma.compress(valid, vari.asma(), axis=-1)

        # Attributes
        set_atts(varo, **atts)

    else: # numpy or ma

        varo = vari[pslice + (valid, )]

    return varo

class _Base_(object):

    def __init__(self, logger=None, **kwargs):

        if logger is not False:

            self.logger = get_logger(logger, **kwfilter(kwargs, 'logger_'))

            self.logger.debug('Instantiate '+self.__class__.__name__)

    def created(self, path):
        """Issue an :meth:`info` message about the creation of a file"""
        path = os.path.abspath(path)
        self.info('Created: '+path)

    def error(self, msg):
        """Raise a :class:`SONARError`"""
        raise SONATError(msg)

    def warning(self, msg):
        """Issue a :class:`SONATWarning`"""
        sonat_warn(msg, stacklevel=3)
    warn = warning

    def notice(self, msg):
        self.logger.notice(msg)

    def info(self, msg):
        self.logger.info(msg)

    def verbose(self, msg):
        self.logger.verbose(msg)

    def debug(self, msg):
        self.logger.debug(msg)

    def _dcache_set_(self, *args, **kwargs):
        if not hasattr(self, '_dcache'):
            self._dcache = {}
        dicttree_set(self._dcache, *args, **kwargs)

    def _dcache_get_(self, *args, **kwargs):
        if not hasattr(self, '_dcache'):
            return
        return dicttree_get(self._dcache, *args, **kwargs)

    def _dcache_clean_(self):
        if hasattr(self, '_dcache'):
            del self._dcache

def rescale_itv(itv, factor=1, target="both", min_width=.1):
    """Scale an interval of floats

    Parameters
    ----------
    itv: tuple of floats + optional bounds
        Interval
    factor: float > 0
        More than one increases the interval
    target: string as one of "both", "min", "max" or "none"
    mindv: float
        Minimal interval absolute width
    """
    if target=="none":
        return itv
    dv = 0.5 * max(itv[1] - itv[0], min_width)
    mean = 0.5 * (itv[0] + itv[1])
    v0 = mean - factor * dv
    v1 = mean + factor * dv
    itvo = [v0, v1]
    if target!="both" and target!="min":
        itvo[0] = itv[0]
    if target!="both" and target!="max":
        itvo[1] = itv[1]
    return tuple(itvo) + itv[2:]


class _XYT_(object):
    """Needs lons, lats and times to be defined"""

    @property
    def ctimes(self):
        if not hasattr(self, '_ctimes'):
            if self.times is None:
                self._ctimes = None
            else:
                self._ctimes = comptime(self.times)
#            self._ctimes.sort()
        return self._ctimes

    def get_lon(self, margin=0):
        """Get the strict (lonmin, lonmax) interval"""
        return rescale_itv((self.lons.min(), self.lons.max()), factor=margin+1)

    def get_lat(self, margin=0):
        """Get the strict (latmin, latmax) interval"""
        return rescale_itv((self.lats.min(), self.lats.max()), factor=margin+1)

    def get_time(self):
        """Get the strict (ctmin, ctmax) interval"""
        return (min(self.ctimes), max(self.ctimes))

    def get_seldict(self, axes='xyt', xybounds='cce', tbounds='cce'):
        sel = {}
        if 'x' in axes:
            sel['lon'] = self.get_lon()
            if xybounds:
                sel['lon'] += xybounds,
        if 'y' in axes:
            sel['lat'] = self.get_lat()
            if xybounds:
                sel['lat'] += xybounds,
        if 't' in axes and self.ctimes:
            sel['time'] = self.get_time()
            if tbounds:
                sel['tile'] += tbounds,
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

class _NamedVariables_(object):

    def check_variables(self, searchmode='ns'):
        """Check that all input variables have known prorties

        See also
        --------
        :func:`sonat.misc.check_variables`
        """
        check_variables([pack.input for pack in self], searchmode=searchmode)

    def select_variables(self, varnames=None, source=None, prefix_to_rm=None):
        """Select variables according to their prefix name

        Parameters
        ----------
        varnames: None, strings
            Selected var names. Defaults to :attr:`varnames`
        source: None, arrays
            Source of array to feed selection. Defaults to :attr:`variables`.
        prefix_to_rm: string
            Prefix to remove before checking id.

        Return
        ------
        list of arrays
        """
        if isinstance(varnames, str):
            varnames = [varnames]
        if source is None:
            source = self.variables
        def isvalid(var):
            id = var.id
            if prefix_to_rm and id.startswith(prefix_to_rm):
                id = id[len(prefix_to_rm):]
            id0 = split_varname(id)[0]
            if varnames is None or id0 in varnames:
                return id
            return False
        vv = []
        for var in source:
            id = isvalid(var)
            if id:
                vv.append(MV2.array(var, id=id, copy=0))
        return vv


    def set_named_norms(self, *anorms, **knorms):
        """Set norms by variable names

        Return
        ------
        list
            The :class:`~sonat.pack.Packer` instance of variables that were
            not normed are returned.

        Example
        -------
        >>> obs.set_named_norms(dict(temp=2.), sal=.8)
        """
        dnorms = dict(*anorms, **knorms)
        notnormed = []
        restack = False
        for packer in self:
            for varname, norm in dnorms.items():
                if (norm is not None and packer.id and
                    packer.id.split('_')[0] == varname):
                    packer.norm = norm
                    restack = True
                    break
            else:
                notnormed.append(packer)
        if restack:
            self._core_stack_()
        return notnormed

def var_prop_dict2list(variables, *args, **kwargs):
    """Convert a variables properties as dict to a list form
    
    If
    
    Parameters
    ----------
    variables: list
        A list of MV2.array or a single array
    \*arnorms, \**kwnorms
        Dictionary arguments whose keys are the :attr:`id` attribute
        prefix of a variable.
        
    Return
    ------
    list
        Properties as a list of the same length as input
        

    Example
    -------
    
    >>> var_prop_dict2list([temp, sal], dict(u=4, temp=6))
    [6, None]
    """
    dargs = dict(*args, **kwargs)
    al = ArgList(variables)
    largs = []
    for var in al.get():
        if hasattr(var, 'id'):
            id = var.id.split('_')[0]
            val = dargs.get(id, None)
        else:
            val = None
        largs.append(val)
    return largs
        

def validate_varnames(varnames):
    """Check that all variable names are in
    :attr:`vacumm.data.cf.GENERIC_VAR_NAMES` list

    Suffixes in the ``_<suffix>`` are removed before searching.

    Raises
    ------
    :class:`sonat.SONATError` if var name is invalid.
    """
    if isinstance(varnames, basestring):
        varnames = [varnames]
    for varname in varnames:
        varname = varname.split('_')[0]
        if varname not in GENERIC_VAR_NAMES:
            raise SONATError('Invalid generic name: '+varname)


def split_varname(varname):
    """Split a variable name in three parts (physical, depth, others)

    The separator is the underscore sign.

    Examples
    --------
    >>> print split_varname('temp_surf_std_dev')
    ('temp', 'surf', 'std_dev')
    >>> print split_varname('temp_variance')
    ('temp', None, 'variance')
    >>> print split_varname('temp')
    ('temp', None, 'None')
    """
    if cdms2.isVariable(varname):
        varname = varname.id
    svname = varname.split('_')
    physical = svname[0]
    depth = others = None
    if len(svname)>1:
        if svname[1] in ['surf', 'bottom']:
            depth = svname[1]
            svname = svname[2:]
            if depth=='3d':
                depth = None
        else:
            svname = svname[1:]
        if svname:
            others = '_'.join(svname)
    return physical, depth, others


def check_variables(vv, searchmode='ns', format=True):
    """Check that all variables are of MV2.array type and are well known
    of the :mod:`vacumm.data.cf` module

    Variable are first checked on the first part of their id, before a "_".
    If the prefix is not known, their are checked is their are known
    with the :func:`vacumm.data.cf.match_known_var`. In this latter case,
    the id of the variable is changed in case of success if format is True.

    Return
    ------
    list
        The corresponding list of generic names
    """
    al = ArgList(vv)
    gennames = []
    for var in al.get():

        # It is MV2.array?
        if not cdms2.isVariable(var):
            raise SONATError('Variable is not of MV2.array type')

        # Is the name generic?
        vns = var.id.split('_')
        varname = vns[0]
        suffix = '_'.join(vns[1:])
        if varname in GENERIC_VAR_NAMES:
            gennames.append(var.id)
            continue

        # Do its properties match a known variable?
        genname = match_known_var(var, searchmode=searchmode)
        if not genname:
            raise SONATError('Unkown variable')
        if format:
            format_var(var, genname, format_axes=False, force=1 if suffix else 2)

        # Suffix
        if suffix:
            if format:
                var.id = var.id + '_' + suffix
            genname = genname + '_' + suffix

        gennames.append(genname)

    return al.put(gennames)

def slice_gridded_var(var, member=None, time=None, depth=None, lat=None, lon=None):
    """Make slices of a variable and squeeze out singletons to reduce it

    The "member" axis is considered here as a generic name for the first
    axis of unkown type.

    .. warning:: All axes must be 1D
    """

    # Check order
    var = var(squeeze=1)
    order = var.getOrder()

    # Unkown axis
    if '-' in order and member is not None:
        i = order.find('-')
        id = var.getAxisIds()[i]
        if isinstance(member, slice):
            kw = {id:member}
            var = var(**kw)
        else:
            axo = create_axis(member)
            cp_atts(var.getAxis(i), axo)
            var = regrid1d(var, axo, iaxi=i)(squeeze=N.isscalar(member))

    # Time interpolation
    if 't' in order and time is not None:
        axi = var.getTime()
        if isinstance(time, slice):
            var = var(time=time)
        else:
            axo = create_time(time, axi.units)
            var = regrid1d(var, axo)(squeeze=N.isscalar(time))

    # Depth interpolation
    if 'z' in order and depth is not None:
        if depth=='bottom':
            var = slice_bottom(var)
        else:
            if depth=='surf':
                depth = slice(-1, None)
            if isinstance(depth, slice):
                var = var(level=depth, squeeze=1) # z squeeze only?
            elif (N.isscalar(depth) and var.getLevel()[:].ndim==1 and
                  depth in var.getLevel()):
                var = var(level=depth)
            else:
                axo = create_dep(depth)
                if axo[:].max()>10:
                    sonat_warn('Interpolation depth is positive. Taking this opposite')
                    axo[:] *=-1
                var = regrid1d(var, axo)(squeeze=N.isscalar(depth))

    # Point
    if (order.endswith('yx') and lon is not None and lat is not None and
            not isinstance(lat, slice) and not isinstance(lon, slice)):

        var = grid2xy(var, lon, lat)(squeeze=N.isscalar(lon))

    else:

        # Latitude interpolation
        if 'y' in order and lat:
            if isinstance(lat, slice):
                var = var(lat=lat)
            else:
                axo = create_lat(lat)
                var = regrid1d(var, axo)(squeeze=N.isscalar(lat))

        # Longitude interpolation
        if 'x' in order and lon:
            if isinstance(lon, slice):
                var = var(lon=lon)
            else:
                axo = create_lon(lon)
                var = regrid1d(var, axo)(squeeze=N.isscalar(lon))

    return var

def slice_bottom(var):
    """Get the bottom value of a variable with 1D z axis

    It takes the deepest unmasked data.
    """
    zz = var.getLevel()
    if zz is None:
        return var
    varo = var(level=slice(0, 1), squeeze=1)
    mask = var.mask
    if not mask.any() or mask.all():
        return varo
    axis = var.getOrder('z')
    for iz in range(1, len(axis)):
        varz = var(level=slice(iz, iz+1), squeeze=1)
        varo[:] = MV2.where(varo.mask, varz, varo, copy=1)
    return varo

def mask_scattered_locs(lons, lats, depths, slice_type, interval, data=None):
    """Maskout coordinates that does not fall within a given interval of lon, lat or depth

    In the case of an vertical interval that varies with location (i.e bottom
    case), no selection is performed but a masking through the data array.

    Parameters
    ----------
    lons: n-D array
    lats: n-D array
    depths: n-D array
    slice_type: one of "zonal", "meridional", "horizontal"
    interval: tuple of float
        Interval for selecting valid data

    Return
    ------
    None or dict
        If not intersect if found, it returns None.
        Else, it returns a dict of lons, lats and depths.

    Todo
    ----
    Add time support.
    """
    # Check slice type
    valid_slice_types = ['zonal', 'merid', 'horiz']
    assert slice_type in valid_slice_types, ('Invalid slice type. '
        'It must be one of: '+', '.join(valid_slice_types))


    # Config
    profiles = (len(depths)!=len(lons) or isaxis(depths) or
                (data is not None and data.ndim==2))
    if cdms2.isVariable(data):
        data = data.asma()
    if data is not None:
        data = data.copy()
        masked_value = True if data.dtype.char=='?' else N.ma.masked

    # Valid locs
    if slice_type is 'horiz':
        bottom = profiles and isinstance(interval[0], N.ndarray)
        if bottom: # varying interval with fixed depth axis
            zz = depths[:, None]
            zz = N.ma.repeat(depths, len(lons), axis=-1) # (nz,np)
        else:
            zz = depths[:] # (nz) or (np)
            if profiles:
                zz = N.ma.repeat(zz.reshape(-1, 1), len(lons), axis=1)
        valid = zz > interval[0]
        valid &= zz <= interval[1]
    else:
        against = lons[:] if slice_type=='merid' else lats[:]
        valid = against >= interval[0]
        valid &= against <= interval[1] # (np)
        if profiles:
            valid = N.resize(valid, depths.shape + valid.shape) # (nz,np)
    mask = ~valid

    # Checks
    if not valid.any():
        return

    # Masking
    if data is not None:
        data[mask] = masked_value
    else:
        data = mask

    return data

def dicttree_relpath(dd, refdir):
    """Make paths stored in a tree of dictionaries relative to another one

    A dict entry must be either a string, a list of strings or another dict
    """
    if isinstance(dd, basestring):
        dd = os.path.relpath(dd, refdir)
    elif isinstance(dd, list):
        for i, d in enumerate(dd):
            dd[i] = os.path.relpath(d, refdir)
    else:
        for key, val in dd.items():
            dd[key] = dicttree_relpath(val, refdir)
    return dd

def interpret_level(level, astuple=False):
    """Interpret a level and format it

    List are returned as lists.
    Non strings are returned whithout change.
    Strings are lower cased.
    Special "surf" and "bottom" values are returned whithout change.
    Others are converted to floats.

    Formats::

        None
        'surf'
        1.3
        [1.2, 5.6]
        (1.3, 'bottom')
        {'var1': ('surf', -53., [-100, -50])}

    """

    # From dict
    if isinstance(level, dict):
        for key, val in level.items():
            if not isinstance(val, tuple):
                val = (val, )
            level[key] = interpret_level(val)
        return level

    # From tuple
    if isinstance(level, tuple):
        out = ()
        for lev in level:
            out += interpret_level(lev),
        return out

    # From list
    if isinstance(level, list):
        for i, d in enumerate(level):
            level[i] = interpret_level(d)

    # Scalar
    elif isinstance(level, basestring):
        level = level.lower()
        if level not in ('bottom', 'surf', '3d'):
            try:
                level = float(level)
            except:
                pass

    return (level, ) if astuple else level


def vminmax(data, asdict=False, symetric=False):
    """Get min and max of dict, list, tuple, array"""
    if isinstance(data, dict):
        vmin, vmax  = vminmax(data.values())
    else:
        vmin, vmax = minmax(data)
    if symetric:
        vmax = min(abs(vmin), abs(vmax))
        vmin = -vmax
    if asdict:
        return dict(vmin=vmin, vmax=vmax)
    return vmin, vmax


def get_long_name(var, default=None):
    """Try to get a long_name"""
    if hasattr(var, 'long_name'):
        return var.long_name
    if hasattr(var, 'id'):
        return var.id.title().replace('_', ' ')
    return default


def recursive_transform_att(data, att, func, *args, **kwargs):
    """Strip ids of a recursive lists of arrays"""
    if not isinstance(data, list):
        if cdms2.isVariable(data) and hasattr(data, att):
                setattr(data, att, func(getattr(data, att), *args, **kwargs))
        return data
    return [recursive_transform_att(dat, att, func, *args, **kwargs) for dat in data]

def sqrt_errors_norm(errors):
    """Compute the norm of square root error (co-)variances

    The output error norm :math:`\sigma` is computed as the square root
    of the mean quadratic errors:

    .. math:: \\sigma = \\sqrt{ \\frac{1}{N} \\sum r^{2}_{i}}

    """
    return N.ma.sqrt((errors**2).mean())
