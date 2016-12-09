"""
Ensembles
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

import re
from collections import OrderedDict
import scipy.stats as SS
from vcmq import (cdms2, MV2, DS, ncget_time, lindates, ArgList,
    MV2_concatenate, create_axis, N, kwfilter, check_case, match_known_var)

from .__init__ import get_logger, sonat_warn, SONATError
from .misc import (list_files_from_pattern, ncfiles_time_indices, asma, NcReader,
    validate_varnames)
from .stack import Stacker
from ._fcore import f_eofcovar, f_sampleens


def load_model_at_regular_dates(ncpat, varnames=None, time=None, lat=None, lon=None,
       level=None,  modeltype='mars', nt=50, dtfile=None, sort=True, asdict=False,
       logger=None, depths=None, **kwargs):
    """Read model output at nearest unique dates with optional linear interpolation


    Parameters
    ----------
    ncpat: string or list of strings
    varnames: string, strings
        Generic var names. If None, all variables that are known from the
        :mod:`vacumm.data.cf` module are used.
    depths: string, list of floats, array, tuple of them, dict

    Examples
    --------
    >>> mdict = load_model_at_regular_dates('myfile.nc', depths='surf')
    >>> mdict = load_model_at_regular_dates('myfile.nc', depths=('surf', 'bottom')
    >>> mdict = load_model_at_regular_dates('myfile.nc', varnames=['temp', 'sal'],
        depths={'temp':('surf', 'bottom'), 'sal':[-50, -10]})

    """
    # Logger
    kwlog = kwfilter(kwargs, 'logger_')
    if logger is None:
        logger = get_logger(**kwlog)
    logger.debug('Loading model at regular dates')


    # Get file list
    ncfiles = list_files_from_pattern(ncpat, time, dtfile=dtfile, sort=True)
    if not ncfiles:
        raise SONATError('No file found')

    # Time interval
    reqtime = time
    if time is None:

        # First
        taxis = ncget_time(ncfiles[0])
        if taxis is None:
            raise SONATError("Can't get time axis for: "+ncfiles[0])
        ctimes = taxis.asComponentTime()
        ct0 = ctimes[0]

        # Last
        if ncfiles[0]!=ncfiles[-1]:
            taxis = ncget_time(ncfiles[-1])
            if taxis is None:
                raise SONATError("Can't get time axis for: "+ncfiles[-1])
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
        raise SONATError(msg)

    # Read
    single = isinstance(varnames, basestring)
    if single:
        varnames = [varnames]
    out = OrderedDict()
    vlevels = {}
    for ncfile, tslices in iidict.items():

        # Dataset instance
        ds = DS(ncfile, modeltype, logger=logger)

        # List of well known variables
        if varnames is None:
            ncvarnames = ds.get_variable_names()
            varnames = []
            for ncvarname in ds.get_variable_names():
                varname = match_known_var(ds[0][ncvarname])
                if varname:
                    varnames.append(varname)

        # Loop on variables
        vardepth = None
        kwvar = dict(lat=lat, lon=lon, verbose=False, bestestimate=False)
        for vname in list(varnames):

            # Level selector
            if vname in vlevels: # cached
                vlevel = vlevels[vname]
            else:
                if isinstance(level, dict): # variable specific
                    vlevel = level.get(vname, None)
                else:
                    vlevel = level # same for all variables
                vlevels[vname] = vlevel # cache it
            if not isinstance(vlevel, tuple):
                vlevel = vlevel,

            # Loop on level specs
            for vlev in vlevel:

                # Output vname
                if isinstance(vlev, basestring) and not vlev=='3d':
                    vnameo = vname + '_' + vlev
                else:
                    vnameo = vname

                # Read and aggregate
                vout = out.setdefault(vnameo, [])
                for tslice in tslices:

                    # Get var
                    kwvar['time'] = tslice
                    var = ds(vname, level=vlev, **kwvar)

                    # Interpolate at depths
                    if ((vlev=='3d' or not isinstance(vlev, basestring)) and
                            depths is not None and var.getLevel() is not None):

                        # Get depths
                        if vardepth is None:
                            vardepth = ds.get_depth(level=vlev, **kwvar)

                        # Interpolate
                        var = ds._interp_at_depths_(var, vardepth, depths)

                    var.id = vnameo
                    vout.append(var)

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
        getmodes=False, logger=None, asdicts=False, anomaly=True, **kwargs):
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
    list (or dict) of arrays:
        variables with their name as keys
    dict: eofs, ev and variance, optional
        eofs: list (or dict) of arrays(nmodes, ...), optional
            EOFs
        ev: array(nmodes), optional
            Eigen values
        var: array
            Variance

    """
    # Logger
    kwlog = kwfilter(kwargs, 'logger_')
    if logger is None:
        logger = get_logger(**kwlog)

    # Ensembe size
    enrich = max(enrich, 1.)
    nt = int(nrens * enrich)

    # Read variables
    data = load_model_at_regular_dates(ncpat, nt=nt, asdict=False, **kwargs)
    single = not isinstance(data, list)

    # Enrichment
    witheofs = nrens!=nt
    if witheofs:

        # Stack packed variables together
        stacker = Stacker(data, norms=norms, logger=logger)
        meanstate = N.zeros(stacker.ns)
        states = N.asfortranarray(stacker.stacked_data.copy())

        # Compute EOFs
        stddev, svals, svecs, status = f_eofcovar(dim_fields=stacker.ns, offsets=1,
            remove_mstate=0, do_mv=0, states=states, meanstate=meanstate)
        if status!=0:
           raise SONATError('Error while calling fortran eofcovar routine')
        neof = svals.size # computed
        neofr = nrens - 1 # retained
        svals = svals[:neofr] * N.sqrt((neof-1.) / neof) # to be consistent with total variance
        svecs = svecs[:, :neofr]

        # Generate ensemble
        sens = f_sampleens(svecs, svals, meanstate, flag=0)

        # Unstack
        ens = stacker.unstack(sens, format=2, rescale='norm' if anomaly else True)
        if getmodes:

            # Modes
            mode_axis = create_axis(N.arange(neofr, dtype='i'), id='mode')
            eofs = stacker.unstack(svecs, firstdims=mode_axis, id='{id}_eof',
                rescale=False, format=1)
            svals = MV2.array(svals, axes=[mode_axis], id='ev',
                attributes={'long_name':'Eigen values'})
            svals.total_variance = float(stacker.ns)

            # Variance
            vv = stacker.format_arrays([d.var(axis=0) for d in stacker.datas],
                id='{id}_variance', mode=1)
            variance = stacker.unmap(vv)

    else: # No enrichment -> take the anomaly if requested

        ens = data

        if anomaly:
            if single:
                ens[:] = ens.asma()-ens.asma().mean(axis=0)
            else:
                for i, e in enumerate(ens):
                    ens[i][:] = e.asma()-e.asma().mean(axis=0)


    # Finalize
    getmodes = getmodes and witheofs
    member_axis = create_axis(N.arange(nrens, dtype='i'), id='member')
    if single:
        ens.setAxis(0, member_axis)
    else:
        for var in ens:
            var.setAxis(0, member_axis)
    if asdicts:
        if single:
            ens = OrderedDict([(ens.id, ens)])
            if getmodes:
                eofs = OrderedDict([(eofs.id, eofs)])
        else:
            ens = OrderedDict([(var.id, var) for var in ens])
            if getmodes:
                eofs = OrderedDict([(var.id, var) for var in eofs])

    # Return
    if not getmodes:
        return ens
    return ens, dict(eofs=eofs, eigenvalues=svals, variance=variance)


class Ensemble(Stacker):
    """Class for exploiting an ensemble"""

    def __init__(self, data, ev=None, variance=None, logger=None, norms=None,
            means=None, **kwargs):

        # Init base
        kwargs['nordim'] = False
        Stacker.__init__(self, data, logger=logger, norms=norms, means=means, **kwargs)
        self.variables = self.inputs

        # Add more variables
        self.ev = ev
        self.variance = variance

    @classmethod
    def from_file(cls, ncfile, varnames=None, ncreader='mars3d',
            exclude=['.*_eof$', 'bounds_.*'],
            include=None, evname='ev', varipat='{varname}_variance',
            lon=None, lat=None, level=None, time=None,
            readertype='generic',
            **kwargs):
        """Init the class with a netcdf file"""

        # Open
        f = NcReader(ncfile, readertype)
        allvars = f.get_variables()

        # Get the list of variables
        varinames = []
        if varnames is None: # guess

            # Basic list
            single = False
            varnames = list(allvars)

            # Filter
            if exclude is None:
                exclude = []
            elif isinstance(exclude, basestring):
                exclude = [exclude]
            if include is None:
                include = []
            elif isinstance(include, basestring):
                    include = [include]
            exclude.append(evname + '$')
            exclude.append('^' + varipat.format(varname='.*') + '$')
            exclude = [re.compile(e, re.I).match for e in exclude]
            include = [re.compile(e, re.I).match for e in include]
            for varname in list(varnames):

                if ((exclude and any([e(varname) for e in exclude])) or
                        (include and not any([e(varname) for e in include]))):

                    varnames.remove(varname) # Remove from state variables

                elif varipat: # Ok, now check variance

                    variname = varipat.format(varname=varname)
                    if variname in allvars:
                        varinames.append(variname)

            # We need at least one
            if not varnames:
                raise SONATError('No valid variable found in file {ncfile}'.format(
                    **locals()))

        else: # provided

            single = isinstance(varnames, N.ndarray)
            if single:
                varnames = [varnames]

            # State variables
            for varname in varnames:
                if varname not in allvars:
                    raise SONATError('Variable {varname} not found in file {ncfile}'.format(
                        **locals()))

            # Variance of state variables
            if varipat:
#                notfound = []
                for varname in varnames:
                    variname = varipat.format(varname=varname)
                    if variname in allvars:
                        varinames.append(variname)
#                    else:
#                        notfound.append(variname)
#                if notfound:
#                    raise SONATError('Variance variables not found: '+' '.join(notfound))

        # Eigen values
        if evname not in allvars:
            evname = None


        # Read
        kwsel = dict(time=time, level=level, lat=lat, lon=lon)
        data = []
        for varname in varnames:
            data.append(f(varname, **kwsel))
        if evname:
            kwargs['ev'] = f('ev', **kwsel)
        if varinames:
            kwargs['variance'] = []
            for variname in varinames:
                vari = f(variname, **kwsel)
                if vari is not None:
                    kwargs['variance'].append(vari)
            kwargs['variance'] = kwargs['variance'] or None
        f.close()

        if single:
            data = data[0]
        return cls(data, **kwargs)

    @property
    def varnames(self):
        return [getattr(input, 'id', None) for input in self.inputs]

    def get_variable(self, varname, depths):
        """Get a variable knowing its name and depths specifications"""
        if hasattr(depths, 'depths'):
            depths = depths.depths
        if isinstance(depths, basestring):
            vname = vname + '_' + dephts
        for var in self.variables:
            if var.id == varname:
                return var
        raise SONATError('Invalid variable name: ' + vname)

    def _ss_diag_(self, diag_name, diag_dict, index=None, format=1, **kwargs):
        """Compute a scipy.stats diagnostic are store it in diags"""
        ss_func = getattr(SS, diag_name)
        if index is not None:
            ss_func_ = ss_func
            ss_func = lambda *args, **kwargs: ss_func_(*args, **kwargs)[index]
        dds = [ss_func(self.datas[i], axis=0, **kwargs) for i in range(self.nvar)]
        dds = self.format_arrays(dds, firstdims=False,
            id='{id}_'+diag_name, mode=format)
        diag_dict[diag_name] = self.unmap(dds)


    def get_diags(self, mean=True, variance=True, kurtosis=True, skew=True,
            skewtest=True, kurtosistest=True, normaltest=True):
        """Get ensemble diagnostics as a dict


        Sea also
        --------
        :func:`scipy.stats.kurtosis` :func:`scipy.stats.kurtosistest`
        :func:`scipy.stats.skew` :func:`scipy.stats.skewtest`
        :func:`scipy.stats.normaltest`
        """
        if (not mean and not variance and not kurtosis and not skewness and
                not skewtest and not kurtosistest and not normaltest):
            return {}
        diags = {}

        # Inputs
        mmeans = self.get_means()
        manoms = [(self.datas[i]-mmeans[i]) for i in range(len(self))]

        # Mean
        if mean is True:
            means = self.fill_arrays(mmeans, firstdims=False,
                id='{id}_mean', format=2)
            diags['mean'] = self.unmap(means)

        # Variance
        if variance is True:

            # Local variance
            vars = [v.var(axis=0) for v in manoms]
            vars = self.fill_arrays(vars, firstdims=False,
                id='{id}_variance', format=2)
            diags['variance'] = self.unmap(vars)

            # Explained variance
            if self.ev is not None and hasattr(self.ev, 'total_variance'):
                diags['explained_variance'] = (self.ev**2).sum()/self.ev.total_variance

            # Local explained variance
            if self.variance is not None:
                lvars = self.remap(self.variance)
                evars = []
                for lvar, var in zip(self.remap(self.variance), vars):
                    if cdms2.isVariable(lvar):
                        evar = lvar.clone()
                        evar.id = var.id.replace('_variance', '_expl_variance')
                    else:
                        evar = lvar.copy()
                    evar[:] /= var



        # Skew
        if skew:
            self._ss_diag_('skew', diags)

        # Kurtosis
        if kurtosis:
            self._ss_diag_('kurtosis', diags)

        # Skew test
        if skewtest:
            self._ss_diag_('skewtest', diags, 0)

        # Kurtosis test
        if kurtosistest:
            self._ss_diag_('kurtosistest', diags, 0)

        # Normal test
        if normaltest:
            self._ss_diag_('normaltest', diags, 0)

        return diags


