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
    MV2_concatenate)

from .__init__ import _Base_, get_logger, pyarm_warn, PyARMError
from .misc import list_needed_files
from .stack import Stacker


def load_model_at_dates(ncpat, ncvars=None, time=None, lat=None, lon=None,
       levels=None,  modeltype='mars', nt=50, dtfile=None, sort=True):
    """Read model output at nearest unique dates with optional linear interpolation

    """

    # Get file list
    ncfiles = list_needed_files(ncpat, time, dtfile=dtfile, sort=True)
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
        if reqdate:
            msg = msg + (", and your requested time range must be enclosed "
                "by model time range.")
        raise PyARMError(msg)

    # Read
    out = OrderedDict()
    vlevels = {}
    for ncfile, slices in iidict.items():

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
            for sl in slices:
                vout.append(ds(vname, time=sl, lat=lat, lon=lon, level=vlevel))

    # Concatenate
    for vname, vlist in vout.items():
        vout[vname] = MV2_concatenate(vlist)

    return vout

def generate_speudo_ensemble(ncpat, nrens=50, enrich=2., norms=None, **kwargs):
    """Generate a static pseudo-ensemble from a single simulation


    Parameters
    ----------
    ncpat: string
        netcdf file name or pattern
    nrens: int
        Ensemble size
    enrich: float
        Enrichment factor
    **kwargs:
        Extra parameters are passed to :func:`load_model_at_dates`

    Return
    ------
    dict: variables with their name as keys

    """

    # Ensembe size
    enrich = max(enrich, 1.)
    nt = int(nrens * enrich)

    # Read variables
    data = load_model_at_dates(ncpat, nt=nt, **kwargs)

    # Enrichment
    if nrens!=nt:

        # Stack packed variables together
        staker = Stacker(data, norms=norms)
        meanstate = N.zeros(stacker.ns)

        # Compute EOFs
        stddev, svals, svecs, status = f_eofcovar(dim_fields=stacker.ns, offsets=1, remove_mstate=0,
            do_mv=0, states=staker.stacked_data, meanstate=meanstate)
        if status!=0:
           raise PyARMError('Error while calling fortran eofcovar routine')

        # Generate ensemble
        sens = f_sampleens(svecs, svals, meanstate, flag)

        # Unstack
        eofs = stacker.unstack(svecs)
        ens = stacker.unstack(sens)
