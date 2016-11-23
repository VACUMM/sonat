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
from vcmq import (cdms2, MV2, DS, ncget_time, lindates,
    MV2_concatenate, create_axis, N, kwfilter)

from .__init__ import _Base_, get_logger, pyarm_warn, PyARMError
from .misc import list_files_from_pattern, ncfiles_time_indices
from .stack import Stacker
from ._fcore import f_eofcovar, f_sampleens


def load_model_at_regular_dates(ncpat, ncvars=None, time=None, lat=None, lon=None,
       level=None,  modeltype='mars', nt=50, dtfile=None, sort=True, asdict=False):
    """Read model output at nearest unique dates with optional linear interpolation

    """

    # Get file list
    ncfiles = list_files_from_pattern(ncpat, time, dtfile=dtfile, sort=True)
    if not ncfiles:
        raise PyARMError('No file found')

    # Time interval
    reqtime = time
    if time is None:

        # First
        taxis = ncget_time(ncfiles[0])
        if taxis is None:
            raise PyARMError("Can't get time axis for: "+ncfiles[0])
        ctimes = taxis.asComponentTime()
        ct0 = ctimes[0]

        # Last
        if ncfiles[0]!=ncfiles[-1]:
            taxis = ncget_time(ncfiles[-1])
            if taxis is None:
                raise PyARMError("Can't get time axis for: "+ncfiles[-1])
            ctimes = taxis.asComponentTime()
        ct1 = ctimes[-1]

        # Time
        time = (ct0, ct1)

    # Generate dates
    dates = lindates(time[0], time[1], nt)

    # Get time indices
    iidict, iiinfo = ncfiles_time_indices(ncfiles, dates, getinfo=True, asslices=True)
    if iiinfo['missed'] or iiinfo['duplicates']:
        msg = ("You must provide at least {nt} model time steps to read "
            "independant dates")
        if reqtime:
            msg = msg + (", and your requested time range must be enclosed "
                "by model time range.")
        raise PyARMError(msg)

    # Read
    single = isinstance(ncvars, basestring)
    if single:
        ncvars = [ncvars]
    out = OrderedDict()
    vlevels = {}
    for ncfile, tslices in iidict.items():

        # Dataset instance
        ds = DS(ncfile, modeltype)

        # List of variables
        if ncvars is None:
            ncvars = [('+'+vname) for vname in ds.get_variable_names()]

        # Loop on variables
        for vname in list(ncvars):

            # Level selector
            if vname in vlevels:
                vlevel = vlevels[vname]
            else:
                if isinstance(level, dict):
                    vlevel = level.get(vname, None)
                else:
                    vlevel = level
                vlevels[vname] = vlevel

            # Read and aggregate
            vout = out.setdefault(vname, [])
            for tslice in tslices:
                vout.append(ds(vname, time=tslice, lat=lat, lon=lon, level=vlevel,
                    verbose=False, bestestimate=False))

    # Concatenate
    for vname, vout in out.items():
        out[vname] = MV2_concatenate(vout)

    # Dict
    if asdict:
        return out

    # Single
    out = out.values()
    if single:
        return out[0]
    return out

def generate_pseudo_ensemble(ncpat, nrens=50, enrich=2., norms=None,
        getmodes=False, logger=None, **kwargs):
    """Generate a static pseudo-ensemble from a single simulation


    Parameters
    ----------
    ncpat: string
        netcdf file name or pattern
    nrens: int
        Ensemble size
    enrich: float
        Enrichment factor
    getmodes: bool
        Get also EOFs end eigen values
    **kwargs:
        Extra parameters are passed to :func:`load_model_at_dates`

    Return
    ------
    dict of arrays:
        variables with their name as keys
    dict of arrays: (nmodes, ...), optional
        EOFs
    array: (nmodes), optiona
        Eigen values

    """
    # Logger
    kwlog = kwfilter(kwargs, 'logger_')
    if logger is None:
        logger = get_logger(**kwlog)

    # Ensembe size
    enrich = max(enrich, 1.)
    nt = int(nrens * enrich)

    # Read variables
    data = load_model_at_regular_dates(ncpat, nt=nt, **kwargs)
    single = not isinstance(data, list)

    # Enrichment
    if nrens!=nt:

        # Stack packed variables together
        stacker = Stacker(data, norms=norms, logger=logger)
        meanstate = N.zeros(stacker.ns)

        # Compute EOFs
        stddev, svals, svecs, status = f_eofcovar(dim_fields=stacker.ns, offsets=1, remove_mstate=0,
            do_mv=0, states=stacker.stacked_data, meanstate=meanstate)
        if status!=0:
           raise PyARMError('Error while calling fortran eofcovar routine')
        svals = svals[:nrens]
        svecs = svecs[:, :nrens]

        # Generate ensemble
        sens = f_sampleens(svecs, svals, meanstate, flag=0)

        # Unstack
        mode_axis = create_axis(N.arange(nrens, dtype='i'), id='mode')
        eofs = stacker.unstack(svecs, firstdims=mode_axis, id='eof_{id}', format=1)
        ens = stacker.unstack(sens, format=2)
        svals = MV2.array(svals, axes=[mode_axis], id='ev',
            attributes={'long_name':'Eigen values'})

        # To dict?
        if not single:
            eofs = OrderedDict([(var.id, var) for var in eofs])
            ens = OrderedDict([(var.id, var) for var in ens])

    # Finalize
    member_axis = create_axis(N.arange(nrens, dtype='i'), id='member')
    if single:
        ens.setAxis(0, member_axis)
    else:
        for var in ens.values():
            var.setAxis(0, member_axis)

    if not nrens!=nt or not getmodes:
        return ens
    return ens, eofs, svals

