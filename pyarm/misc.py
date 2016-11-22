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
import codecs
import string
from string import Formatter
from glob import has_magic, glob
from collections import OrderedDict
from vcmq import (ncget_time, itv_intersect, pat2freq, lindates, adatetime,
    comptime, add_time, pat2glob, are_same_units)

from .__init__ import pyarm_warn, PyARMError

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
        - with_time: list of keys that date pattern format
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

                pyarm_warn("You should better provide a time interval "
                    "with a date pattern in file name")
                ncfile = DatePat2GlobFormatter().format(ncpat, **subst)
                files = glob(ncfile)

            else: # with

                # Guess pattern and frequency
                date_format = scan_fields[scan_props['with_time'][0]]['format_spec']
                freq = pat2freq(date_format)
                if dtfile is None:
                    dtfile = 1, freq
                    pyarm_warn('Time steps between files not explicitly specified. '
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

def ncfiles_time_indices(ncfiles, dates, getinfo=False):
    """Get time indices corresponding to each dates for each files

    The time axis is read only for the first two files.
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
                PyARMError("Can't read time axis in file: " + ncfile)
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

    if not getinfo:
        return ncfdict
    return ncfdict, {'missed':dates, 'duplicates':duplicate_dates}
