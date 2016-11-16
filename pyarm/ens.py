"""
"""
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
    iidict, iiinfo = ncfiles_time_indices(ncfiles, dates, getinfo=True, asslices=True))
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
            do_mv=0, staker.stacked_data, meanstate=meanstate)
        if status!=0:
           raise PyARMError('Error while calling fortran eofcovar routine')

        # Generate ensemble
        sens = f_sampleens(svecs, svals, meanstate, flag)

        # Unstack
        eofs = stacker.unstack(svecs)
        ens = stacker.unstack(sens)
    