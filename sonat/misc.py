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

import re
import os
import codecs
import string
from string import Formatter
from glob import has_magic, glob
from collections import OrderedDict
import numpy as N
import cdms2
from vcmq import (ncget_time, itv_intersect, pat2freq, lindates, adatetime,
    comptime, add_time, pat2glob, are_same_units, indices2slices,
    kwfilter, numod, GENERIC_VAR_NAMES, DS, set_atts, format_var,
    match_known_var, ArgList, create_lon, regrid1d, grid2xy, create_lat,
    create_time, create_dep, create_axis, cp_atts)

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
    """List files possibly with glob and date patterns"""

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
    files = [ncfile for ncfile in files if os.path.exists(ncfile)]

    # Unique
    files = list(set(files))

    # Sort
    if sort:
        files.sort(key=sort if callable(sort) else None)

    return files

def ncfiles_time_indices(ncfiles, dates, getinfo=False, asslices=False):
    """Get time indices corresponding to each dates for each files

    .. warning:: The time axis is read only for the first two files.
        All files must be obey the following rules:

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
            ctimes = taxis.asComponentTime()

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
    single = len(args)==1
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

    >>> f = NcReader(ncfile, 'netcdf'') # python netcdf4 reader
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
    """Keep valid spatial points"""

    # Slice for extra dimensions
    pslice = (slice(None), ) * (valid.ndim-2)

    if cdms2.isVariable(vari):

        nv = valid.sum()


        if vari.getGrid() is not None and len(vari.getGrid().shape)==2:

            # Init
            assert valid.ndim == 2, 'Valid must be a 2D array'
            varo = vari[pslice + (0, slice(0, nv))].clone()
            ax = create_axis((nv, ), id='point', long_name='Spatial points')
            varo.setAxis(-1, ax)
            varo.setGrid(None)

            # Fill
            varo[:] = vari.asma()[pslice, valid]

        else:

            # Init
            ax = vari.getAxis(-1)
            varo = vari[pslice + (slice(0, nv), )].clone()

            # Fill
            varo.getAxis(-1)[:] = N.compress(valid, ax[:])
            varo[:] = N.compress(valid, vari.asma(), axis=-1)

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

    def warn(self, msg):
        """Issue a :class:`SONATWarning`"""
        sonat_warn(msg, stacklevel=3)

    def error(self, msg):
        """Raise a :class:`SONARError`"""
        raise SONARError(msg)

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

    def get_seldict(self, axes='xyt', xybounds='cce', tbounds='cce'):
        sel = {}
        if 'x' in axes:
            sel['lon'] = (self.lons.min(), self.lons.max(), xybounds)
        if 'y' in axes:
            sel['lat'] = (self.lats.min(), self.lats.max(), xybounds)
        if 't' in axes and self.ctimes:
            sel['time'] = (min(self.ctimes), max(self.ctimes), tbounds)
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

class _CheckVariables_(object):

    def check_variables(self, searchmode='ns'):
        """Check that all input variables have known prorties

        See also
        --------
        :func:`sonat.misc.check_variables`
        """
        check_variables([pack.input for pack in self], searchmode=searchmode)


def validate_varnames(varnames):
    """Check that all variable names are in
    :attr:`vacumm.data.cf.GENERIC_VAR_NAMES` list

    Suffixes in the ``_<suffix>`` are removed before searching.

    Raise
    -----
    :class:`sonat.SONATError` if var name is invalid.
    """
    if isinstance(varnames, basestring):
        varnames = [varnames]
    for varname in varnames:
        varname = varname.split('_')[0]
        if varname not in GENERIC_VAR_NAMES:
            raise SONATError('Invalid generic name: '+varname)

def check_variables(vv, searchmode='ns', format=True):
    """Check that all variables are of MV2.array type and is well known
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
        if isinstance(depth, slice):
            var = var(level=depth, squeeze=1)
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

#def _take_first_index_(var,
