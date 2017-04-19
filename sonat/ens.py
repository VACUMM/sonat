# -*- coding: utf8 -*-
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

import os
import re
from collections import OrderedDict
import scipy.stats as SS
from matplotlib.ticker import MultipleLocator
from vcmq import (cdms2, MV2, DS, ncget_time, lindates, ArgList, format_var,
    MV2_concatenate, create_axis, N, kwfilter, check_case, match_known_var,
    CaseChecker, curve, bar, dict_check_defaults, dict_merge, checkdir,
    dicttree_get, dicttree_set, latlab, lonlab, P, regrid2d, Plot2D)

from .__init__ import get_logger, sonat_warn, SONATError
from .misc import (list_files_from_pattern, ncfiles_time_indices, asma, NcReader,
    validate_varnames, _NamedVariables_, check_variables, dicttree_relpath,
    interpret_level, split_varname)
from .stack import Stacker
from ._fcore import f_eofcovar, f_sampleens
from .plot import (plot_gridded_var, DEFAULT_PLOT_KWARGS)
from .render import register_html_template, render_and_export_html_template

#: Default plotting keyword per ensemble diagnostic
DEFAULT_PLOT_KWARGS_PER_ENS_DIAG = {
    'spectrum':dict(
        xmin=.5, ymin=0,
        title='Explained variance', width=.5,
        xlabel="%(xlong_name)s",
        bottom=.2,
    ),
    'explained_variance':{
        'cmap':'positive',
        'levels_mode':'positive',
    },
    'variance':{
        'cmap':'positive',
        'levels_mode':'positive',
    },
    'skew':{
        'cmap':'symetric',
        'levels_mode':'symetric',
    },
    'kurtosis':{
        'cmap':'symetric',
        'levels_mode':'symetric',
    },
    'skewtest':{
        'cmap':'symetric',
        'levels_mode':'symetric',
    },
    'kurtosistest':{
        'cmap':'symetric',
        'levels_mode':'symetric',
    },
    'normaltest':{
        'cmap':'positive',
        'levels_mode':'positive',
    },
}



def load_model_at_regular_dates(ncpat, varnames=None, time=None, lat=None, lon=None,
       level=None, depths=None,  modeltype='mars', nt=50, dtfile=None, sort=True, asdict=False,
       logger=None, **kwargs):
    """Read model output at nearest unique dates with optional linear interpolation


    Parameters
    ----------
    ncpat: string or list of strings
    varnames: string, strings
        Generic var names. If None, all variables that are known from the
        :mod:`vacumm.data.cf` module are used.
    level: string, None, list of floats, array, tuple of them, dict
        Here are some possible values:

        - "surf" or "bottom": self explenatory
        - None or "3d": No slice, so get all levels with no interpolation.
        - A list or array of negative depths: get all levels and
          interpolate at these depths.

        Variables sliced with "surf" and "bottom" are returned with
        an id suffixed with "_surf" or "_bottom".
        You can speficy different slicings  using a tuple
        of depth specifications.
        You can specialise slicings of a variable using a dictionary with
        the key as the variable name.

    See also
    --------
    :func:`sonat.misc.list_files_from_pattern` for more options


    Examples
    --------
    >>> mdict = load_model_at_regular_dates('myfile.nc', level='surf')
    >>> mdict = load_model_at_regular_dates('myfile.nc', level=('surf', 'bottom')
    >>> mdict = load_model_at_regular_dates('myfile.nc', varnames=['temp', 'sal'],
        level={'temp':('surf', 'bottom'), 'sal':[-50, -10]})
    >>> mdict = load_model_at_regular_dates('myfile.nc', varnames=['temp', 'sal'],
        level={'temp':('surf', '3d'), 'sal':None}, depths=[-50, -10])


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
    if not isinstance(level, dict):
        level = {'__default__':level}
    for ncfile, tslices in iidict.items():

        # Dataset instance
        ds = DS(ncfile, modeltype, logger_name='SONAT.Dataset', logger_level='error')

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

            # Level selector for this variable
            if vname in vlevels: # cached
                vlevel = vlevels[vname]
            else:
                vlevel = interpret_level(dicttree_get(level, vname), astuple=True)
                vlevels[vname] = vlevel # cache it

            # Loop on level specs
            for vlev in vlevel:

                # Output vname and vlev check
                if not isinstance(vlev, basestring):
                     vnameo = vname
                elif vlev not in ('surf', "bottom", "3d"):
                    raise SONATError('Depth string must one of '
                        'surf, bottom, 3d')
                elif vlev!='3d':
                    vnameo = vname + '_' + vlev
                else:
                    vlev = None
                    vnameo = vname

                # Slicing level and output depths
                if vlev not in ['surf', 'bottom']:

                    # numeric so interpolation
                    if vlev is None:
                        vdep = depths if depths is not None else None
                    else:
                        vdep = vlev
                    interp = vdep is not None
                    if interp:
                        vlev = None

                else:

                    interp = False

                # Read and aggregate
                vout = out.setdefault(vnameo, [])
                vinterp = None
                for tslice in tslices:

                    # Get var
                    kwvar['time'] = tslice
                    var = ds(vname, level=vlev, **kwvar)

                    # Interpolate at numeric depths
                    if interp and var.getLevel() is not None:

                        # Get depths
                        if True or vardepth is None: #FIXME: bad to always read it
                            vardepth = ds.get_depth(level=vlev, zerolid=True,
                                **kwvar)

                        # Interpolate
                        var = ds._interp_at_depths_(var, vardepth, vdep,
                            extrap='top')

                    # Id with suffix
                    var.id = vnameo

                    # Store results
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

def generate_pseudo_ensemble(ncpat, varnames=None, nrens=50, enrich=2., norms=None,
        getmodes=False, logger=None, asdicts=False, anomaly=True, ncensfile=None,
        **kwargs):
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
    logger.verbose('Generating pseudo-ensemble')

    # Ensembe size
    enrich = max(enrich, 1.)
    nt = int(nrens * enrich)
    logger.debug(' enrich={enrich},  nt={nt}, ncpat={ncpat}, varnames={varnames}'.format(**locals()))

    # Read variables
    logger.debug('Reading the model at {} dates'.format(nt))
    data = load_model_at_regular_dates(ncpat, varnames=varnames, nt=nt,
        asdict=False, **kwargs)
    single = not isinstance(data, list)

    # Enrichment
    witheofs = nrens!=nt
    if witheofs:
        logger.debug('Computing reduced rank ensemble with EOFs analysis')

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
            mode_axis = create_axis(N.arange(1, neofr+1, dtype='i'), id='mode')
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

        logger.debug('Getting the anomaly to build the ensemble')
        ens = data

        if anomaly:
            if single:
                ens[:] = ens.asma()-ens.asma().mean(axis=0)
            else:
                for i, e in enumerate(ens):
                    ens[i][:] = e.asma()-e.asma().mean(axis=0)


    # Finalize
    getmodes = getmodes and witheofs
    member_axis = create_axis(N.arange(nrens, dtype='i'), id='member',
        long_name='Member')
    if single:
        ens.setAxis(0, member_axis)
    else:
        for var in ens:
            var.setAxis(0, member_axis)

    # Dump to file
    if ncensfile:
        logger.debug('Dump the ensemble to netcdf')
        checkdir(ncensfile)
        f = cdms2.open(ncensfile, 'w')
        ensvars = list(ens) if not single else [ens]
        if getmodes:
            if single:
                ensvars.append(eofs)
                ensvars.append(variance)
            else:
                ensvars.extend(eofs)
                ensvars.extend(variance)
            ensvars.append(svals)
        for var in ensvars:
            f.write(var)
        f.close()
        logger.created(ncensfile)

    # As dicts
    if asdicts:
        if single:
            ens = OrderedDict([(ens.id, ens)])
            if getmodes:
                eofs = OrderedDict([(eofs.id, eofs)])
                variance = OrderedDict([(variance.id, variance)])
        else:
            ens = OrderedDict([(var.id, var) for var in ens])
            if getmodes:
                eofs = OrderedDict([(var.id, var) for var in eofs])
                variance = OrderedDict([(var.id, var) for var in variance])


    # Return
    if not getmodes:
        return ens
    return ens, dict(eofs=eofs, eigenvalues=svals, variance=variance)


class Ensemble(Stacker, _NamedVariables_):
    """Class for exploiting an ensemble


    Parameters
    ----------
    data: array or list of arrays
        Input model variable with same first axis
    ev: 1D array, optional
        Eigen values from EOF decomposition
    variance: array or list of arrays
        Variance of each variable
    norms: list or dict of floats
        Coefficients used to normalise variance.
        Either a list of ordered coefficient, one per variable,
        or a dict whose keys are variable prefix name.
        .. note:: Norms are set after syncing if sync_norms is True.
    sync_norms: bool, optional
        Make sur that all variable whose id has the same prefix (string
        before "_") have the same norm
    """

    def __init__(self, data, ev=None, variance=None, logger=None, norms=None,
            means=None, sync_norms=True, check_vars=False, bathy=None, **kwargs):

        # Check input variables
        datas = ArgList(data).get()
        if check_vars:
            check_variables(datas)
        for var in datas: # rename first axis
            var.getAxis(0).id = 'member'

        # Init base
        kwargs['nordim'] = False
        Stacker.__init__(self, data, logger=logger,
                         norms=None, means=means, **kwargs)
        self.variables = self.inputs

        # Synchronise norms
        self._norm_synced = False
        if sync_norms:
            self.sync_norms()

        # Norms
        if norms:
            if isinstance(norms, dict):
                self.set_named_norms(norms)
            else:
                self.norms = norms

        # Add more variables
        self.ev = ev
        self.variance = variance

        # Caches
        self._diags = {}
        self._meananoms = {}

        # Bathy
        if bathy is not None:
            self.set_bathy(bathy)


    @classmethod
    def from_file(cls, ncfile, varnames=None, ncreader='mars3d',
            exclude=['.*_eof$', 'bounds_.*'],
            include=None, evname='ev', varipat='{varname}_variance',
            lon=None, lat=None, level=None, time=None,
            readertype='generic', checkvars=False,
            **kwargs):
        """Init the class with a netcdf file"""

        # Open
        f = NcReader(ncfile, readertype, logger_level='error')
        try:
            allvars = f.get_variables()
        except:
            raise SONATError('Error while retreiving the list of variable in '
                             'netcdf file: '+ncfile)

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

        # Check that variables are well known
        if checkvars:
            gennames = check_variables([f[varname] for varname in varnames], format=False)
        kwargs['checkvars'] = False

        # Eigen values
        if evname not in allvars:
            evname = None


        # Read
        kwsel = dict(time=time, level=level, lat=lat, lon=lon)
        data = []
        for i, varname in enumerate(varnames):
            var = f(varname, **kwsel)
            data.append(var)
            if checkvars:
                vns = gennames[i].split('_')
                suffix = '_'.join(vns[1:])
                format_var(var, vns[0], format_axes=False, force=True)
                if suffix:
                    var.id = var.id + '_' + suffix
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
        ens = cls(data, **kwargs)
        ens.debug('From file: '+ncfile)
        ens.ncfile = ncfile
        return ens

    @property
    def varnames(self):
        return [getattr(input, 'id', None) for input in self.inputs]

    def get_prefixed(self, prefix):
        """Get index of all variables that start with prefix"""
        prefix = prefix.strip('_')
        return [i for i, varname in enumerate(self.varnames)
            if varname and varname.split("_")[0].startswith(prefix)]

    def get_prefixes(self):
        return set([varname.split("_")[0] for varname in self.varnames
            if varname])


    def get_variable(self, vname, depths):
        """Get a variable knowing its name and the depths specifications"""
        if hasattr(depths, 'depths'):
            depths = depths.depths
        if isinstance(depths, basestring):
            vname = vname + '_' + depths
        for var in self.variables:
            if var.id == vname:
                return var
        raise SONATError('Invalid variable name: ' + vname)

    def get_grid(self):
        return self.variables[0].getGrid()

    def set_bathy(self, bathy2d):
        """Interpolate a bathymetry and save it"""
        if bathy2d is not None:
            self._bathy = regrid2d(bathy2d, self.get_grid())

    def get_bathy(self, bathy2d=None, force=False):
        if bathy2d is not None and (force or not hasattr(self, '_bathy')):
            self.set_bathy(bathy2d)
        return getattr(self, '_bathy', None)

    bathy = property(fget=get_bathy, fset=set_bathy, doc='Colocated bathymetry')


    def sync_norms(self, force=True):
        """Synchronise norms between variables whose id has the same prefix

        Parameters
        ----------
        force: bool
            Force sync even if it seems at has already been synced.
        """
        if not force and self._norm_synced:
            return False

        # Renorm for each prefix
        for prefix in self.get_prefixes():

            # Indices
            indices = self.get_prefixed(prefix)

            # Unified norm
            psize = sum([self[i].psize for i in indices])
            norm = N.sqrt(sum([(self[i].psize * self[i].norm**2)
                for i in indices]) / psize)

            # Set it
            for i in indices:
                self[i].set_norm(norm)

        self._norm_synced = True
        return True

    def get_named_norms(self, sync=None):
        """Get the norm of each variable prefix

        Example
        -------
        >>> dnorms = ens.get_named_norms(sync=True)
        >>> temp_norm = dnorms['temp']
        """
        # Norms must all be sync by variable prefix
        if sync is None:
            sync = not  self._norm_synced
        if sync:
            self.sync_norms()

        # Fill norms
        dnorms = {}
        for prefix in self.get_prefixes():

            # Indices
            indices = self.get_prefixed(prefix)

            # Get the norm
            dnorms[prefix] = self[indices[0]].norm

        return dnorms


    def _ss_diag_(self, diag_name, diag_dict, index=None, format=1, **kwargs):
        """Compute a scipy.stats diagnostic are store it in diags"""
        # From cache
        if diag_name in self._diags:
            diag_dict[diag_name] = self._diags[diag_name]
            return

        # Compute
        ss_func = getattr(SS, diag_name)
        if index is not None:
            ss_func_ = ss_func
            ss_func = lambda *args, **kwargs: ss_func_(*args, **kwargs)[index]
        dds = [ss_func(self.datas[i], axis=0, **kwargs) for i in range(self.nvar)]
        dds = self.format_arrays(dds, firstdims=False,
            id='{id}_'+diag_name, mode=format)
        for dd in dds:
            dd.sonat_ens_diag = diag_name
        diag_dict[diag_name] = self.unmap(dds)


    def get_diags(self, mean=True, variance=True, kurtosis=True, skew=True,
            skewtest=True, kurtosistest=True, normaltest=True):
        """Get ensemble diagnostics as a dict

        Example
        -------
        >>> diags = ens.get_diags(mean=False, skew=True)
        >>> assert diags['skew'].shape == diags['skewtest'].shape

        See also
        --------
        :func:`scipy.stats.skew`
        :func:`scipy.stats.skewtest`
        :func:`scipy.stats.kurtosis`
        :func:`scipy.stats.kurtosistest`
        :func:`scipy.stats.normaltest`
        """
        if (not mean and not variance and not kurtosis and not skewness and
                not skewtest and not kurtosistest and not normaltest):
            return {}
        diags = OrderedDict()
        self.verbose('Getting ensemble diagnostics')

        # Inputs
        if not self._meananoms:
            mmeans = self.get_means()
            manoms = [(self.datas[i]-mmeans[i]) for i in range(len(self))]
            self._meananoms = {'means':mmeans, 'anoms':manoms}
        else:
            mmeans = self._meananoms['means']
            manoms = self._meananoms['anoms']

        # Mean
        if mean is True:
            if 'mean' in self._diags:
                diags['mean'] = self._diags['mean']
            else:
                means = self.fill_arrays(mmeans, firstdims=False,
                    id='{id}_mean', format=2)
                diags['mean'] = self._diags['mean'] = self.unmap(means)

        # Variance
        if variance is True:

            # Local variance
            if 'variance' in self._diags:
                diags['variance'] = self._diags['variance']
            else:
                vars = [v.var(axis=0) for v in manoms]
                vars = self.fill_arrays(vars, firstdims=False,
                    id='{id}_variance', format=2)
                for var in vars:
                    var.sonat_ens_diag = 'variance'
                diags['variance'] = self._diags['variance'] = self.unmap(vars)



            # Spectrum
            if self.ev is not None and hasattr(self.ev, 'total_variance'):
                if 'spectrum' in self._diags:
                    diags['spectrum'] = self._diags['spectrum']
                else:
                    diags['spectrum'] = self._diags['spectrum'] = \
                        100 * (self.ev**2) / self.ev.total_variance
                    diags['spectrum'].id = 'spectrum'
                    diags['spectrum'].units = '%'
                    diags['spectrum'].sonat_ens_diag = 'spectrum'
                    diags['spectrum'].long_name = 'Relative spectrum'

            # Local explained variance
            if self.variance is not None:
                if 'explained_variance' in self._diags:
                    diags['explained_variance'] = self._diags['explained_variance']
                else:
                    evars = []
                    for lvar, var in zip(diags['variance'], self.remap(self.variance)):
                        if cdms2.isVariable(lvar):
                            evar = lvar.clone()
                            evar.id = var.id[:-len('_variance')] + '_expvar'
                            evar.units = '%'
                        else:
                            evar = lvar.copy()
                        try:
                            evar[:] /= var
                        except:
                            pass
                        evar[:] *= 100
                        evar.sonat_ens_diag = 'explained_variance'
                        evars.append(evar)
                    diags['explained_variance'] = self._diags['explained_variance'] = \
                        self.unmap(evars)


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

    def plot_diags(self, mean=True, variance=True, kurtosis=True, skew=True,
            skewtest=True, kurtosistest=True, normaltest=True,
            titlepat = '{var_name} - {diag_long_name} - {slice_loc}',
            surf=None, bottom=None, horiz_sections=None, #points=None,
            zonal_sections=None, merid_sections=None,
            fmtlonlat='{:.2f}', fmtdep='{:04.0f}m',
            figpat_slice='sonat.ens.{diag_name}_{var_name}_{slice_type}_{slice_loc}.png',
            figpat_generic='sonat.ens.{diag_name}.png',
            show=False, savefig=True,
            props=None, **kwargs):
        """Create figures for diagnostics

        Diagnostic arguments are passed to :meth:`get_diags`,
        like ``variance=True``.

        Parameters
        ----------
        surf: bool
        bottom: bool
        horiz_sections: floats
            List of depths for horizontal slices. All by default.
            It accepts scalars or list of floats and of the two special values
            "bottom" and "surf". In the latter case, these variable must
            be explicitely included in the ensemble, like "temp_surf". A depth
            of 0 is not equivalent to "surf", despite results may be similar,
            or equivalent when the sea level is always 0.
            In the case of floats, the
            ensemble must contain 3d variables.
        zonal_sections: None, list, floats
            Latitudes
        merid_sections: None, list, floats
            Longitudes
        props: dict, None
            Dictionary of graphic properties passed as keywords to plotting
            functions. The keys must be valid diagnostics names such as
            "mean", or one of the follow plotting function: map, curve.

        Return
        ------
        dict
            A tree of :class:`dict` with the following branch types,
            with all keys capitalised:

            - Diagnostic
            - Variable
            - Type of slice
            - Location of slice
        """

        # Diags
        diags = self.get_diags(mean=mean, variance=variance, kurtosis=kurtosis, skew=skew,
            skewtest=skewtest, kurtosistest=kurtosistest, normaltest=normaltest)

        # Inits
        self.verbose('Plotting ensemble diagnostics')
        figs = OrderedDict()
        if props is None:
            props = {}
        kwprops = dict_merge(props, DEFAULT_PLOT_KWARGS_PER_ENS_DIAG)

        # Loop on diags
        figs = OrderedDict()
        for diag_name, diag_var in diags.items():

            diag_long_name = diag_name.title().replace('_', ' ')
            self.debug('Plotting ens diag: '+diag_long_name)

            # Explained variance
            if diag_name=='spectrum':

                self.debug(' Variable: spectrum')
                figfile = figpat_generic.format(**locals())
                checkdir(figfile)
                kw = kwprops.get(diag_name, {})
                dict_check_defaults(kw, xmin=.5, xmax=len(diag_var)+.5,
                                    width=.5, xminor_locator=False,
                                    ymax=1.1*diag_var.max(),
                                    savefig=figfile,
                                    fig='new', xlabel='%(xlong_name)s',
                                    xlocator=MultipleLocator(1),
                                    **DEFAULT_PLOT_KWARGS)
                p = bar(diag_var, **kw)
                self.created(figfile)
                figs[diag_long_name] = figfile

            # Plot other sliced variables
            else:

                ff = self.plot_fields(diag_var,
                                      surf=surf, bottom=bottom,
                                      horiz_sections=horiz_sections,
                                      zonal_sections=zonal_sections,
                                      merid_sections=merid_sections,
                                      points=points,
                                      figpat=figpat_slice, savefig=savefig,
                                      titlepat=titlepat,
                                      fmtlonlat=fmtlonlat, fmtdep=fmtdep,
                                      subst={'diag_long_name':diag_long_name,
                                             'diag_name':diag_name},
                                      props=props,
                                      )
                if ff:
                    figs[diag_long_name] = ff


        return figs


    def plot_fields(self, variables,
                    surf=None, bottom=None, horiz_sections=False,
                    points=None, zonal_sections=None,
                    merid_sections=None, titlepat="{var_name} - {slice_loc}",
                    props=None, subst={},
                    figpat='sonat.ens.{var_name}_{slice_type}_{slice_loc}.png',
                    fmtlonlat=u'{:.2f}°{}', fmtdep='{:04.0f}m',
                    savefig=True, obs=None,
                    close=True, show=False,
                    **kwargs):
        """Plot one or several platform like variables

        Parameters
        ----------
        depths: string, floats
            List of depths for horizontal slices. All by default.
            It accepts scalars or list of floats and of the two special values
            "bottom" and "surf". In the latter case, these variable must
            be explicitely included in the ensemble, like "temp_surf". A depth
            of 0 is not equivalent to "surf", despite results may be similar,
            or equivalent when the sea level is always 0.
            In the case of floats, the
            ensemble must contain 3d variables.
        props: dict, None
            Dictionary of graphic properties passed as keywords to plotting
            functions. The keys must be valid diagnostics names such as
            "mean", or one of the follow plotting function: map, curve.

        """
        # Init
        if props is None:
            props = {}
        kwprops = dict_merge(props, DEFAULT_PLOT_KWARGS_PER_ENS_DIAG)
        kwprops.update(fig='new') # TODO: must be more flexible like for obs
        figs = OrderedDict()
        kwobs = kwfilter(kwargs, 'obs')
        kwobs.update(surf=False, bottom=False, zonal_sections=None,
                     merid_sections=None, horiz_sections=None)

        # Loop on variables
        for variable in ArgList(variables).get():
            var_short_name, var_depth, diag_name = split_varname(variable)
            diag_name = getattr(variable, 'sonat_ens_diag', diag_name)
            var_name = (var_short_name if var_depth is None else
                        (var_short_name + '_' + var_depth))
            varname = var_name
            id = variable.id
            self.debug(' Variable: '+id)
            order = variable.getOrder()
            if diag_name:
                kw = kwprops.get(diag_name, {})
            else:
                kw = {}
            toplot = []

            # Maps
            if (horiz_sections is not False or surf is not False or
                bottom is not False):

                # Surf/bottom/3D
                if ((var_depth=='surf' and surf is not False) or
                        (var_depth=='bottom' and bottom is not False)): # 2D
                    slice_loc = var_depth
                    toplot.append(dict(slice_loc=slice_loc, kw=kw,
                                       keys=(var_short_name, 'map', slice_loc),
                                       obs_plot={var_depth:True},
                                       slice_type = 'map'
                                       ))

                elif (variable.getLevel() is not None and
                      horiz_sections is not False and horiz_sections is not None): # 3D


                    # Surf and bottom
                    for ok, slice_loc in [(surf, 'surf'), (bottom, 'bottom')]:
                        if ok is True:
                            toplot.append(dict(
                                slice_loc=slice_loc,
                                slice_type='map',
                                keys=(var_short_name, 'map', slice_loc),
                                kw=dict_merge(dict(depth=slice_loc), kw),
                                obs_plot={slice_loc:True},
                                ))


                    # Loop on 3D depths specs
                    if horiz_sections is True:
                        depths =  variable.getLevel()
                    else:
                        depths = horiz_sections
                        if N.isscalar(depths):
                            depths = [depths]
                    depths = list(depths)
                    for depth in depths:
                        if isinstance(depth, str):
                            slice_loc = depth
                        else:
                            slice_loc = fmtdep.format(abs(depth))
                        toplot.append(dict(
                            slice_loc=slice_loc,
                            slice_type='map',
                            keys=(var_name, 'map', slice_loc),
                            kw=dict_merge(dict(depth=depth), kw),
                            obs_plot={'horiz_sections':depth},
                            ))

            # Zonal sections
            if zonal_sections is not None and zonal_sections is not False:
                slice_type = 'zonal'
                if N.isscalar(zonal_sections):
                    zonal_sections = [zonal_sections]
                for lat in zonal_sections:
                    slice_loc = latlab(lat, no_symbol=True)
                    toplot.append(dict(
                        slice_loc=slice_loc,
                        slice_type=slice_type,
                        keys=(var_name, slice_type, slice_loc),
                        kw=dict_merge(dict(lat=lat), kw),
                        obs_plot={'zonal_sections':lat},
                        ))

            # Meridional sections
            if merid_sections is not None and merid_sections is not False:
                slice_type = 'merid'
                if N.isscalar(merid_sections):
                    merid_sections = [merid_sections]
                for lon in merid_sections:
                    slice_loc = lonlab(lon, no_symbol=True)
                    toplot.append(dict(
                        slice_loc=slice_loc,
                        slice_type=slice_type,
                        keys=(var_name, slice_type, slice_loc),
                        kw=dict_merge(dict(lon=lon), kw),
                        obs_plot={'merid_sections':lon},
                        ))

#            # Points
#            if points:
#                slice_type = 'point'
#                if isinstance(points, tuple):
#                    points = [points]
#                for lon, lat in points:
#                    slice_loc = (lonlab(lon, no_symbol=True) + '-' +
#                                 latlab(lat, no_symbol=True))
#                    toplot.append(dict(
#                        slice_loc=slice_loc,
#                        slice_type=slice_type,
#                        keys=(var_name, slice_type, slice_loc),
#                        kw=dict_merge(dict(lon=lon, lat=lat), kw),
#                        obs_plot={'points':(lon, lat)},
#                        ))


            # Make all pending plots for this variable
            for spec in toplot:

                # Get specs
                slice_loc = spec['slice_loc']
                slice_type = spec['slice_type']
                obs_plot = spec['obs_plot']
                keys = map(str.title, spec['keys'])
                kw = kwargs.copy()
                kw.update(spec['kw'])
                self.debug('  Location: '+slice_loc)

                # Setup
                dfmt = locals()
                dfmt.update(subst)
                title = titlepat.format(**dfmt)
                kw.update(title=title, savefig=False, close=False)
                dict_check_defaults(kw, fig='new',
                                    **DEFAULT_PLOT_KWARGS_PER_ENS_DIAG)

                # Plot
                p = plot_gridded_var(variable, **kw) # Action !

                # Plot obs locations on 2D plots only
                if obs and isinstance(p, Plot2D):
                    #FIXME: should we pass lon/lat/level from variable?
                    kwo = kwobs.copy()
                    kwo.update(obs_plot, full3d=False, full2d=False, plotter=p,
                               savefig=False, close=False, title=False,
                               colorbar=False)
#                    obs.set_cached_plot('locations', slice_type, slice_loc, p)
                    obs.plot('locations', **kwo)

                # Finalise
                if savefig:
                    figfile = figpat.format(**dfmt)
                    checkdir(figfile)
                    p.savefig(figfile)
                    close = True
                    self.created(figfile)
                    kw = {keys[-1]:figfile, '__class__':OrderedDict}
                    dicttree_set(figs, *keys[:-1], **kw)
                if not show and close:
                    p.close()
        if show:
            P.show()
        elif close:
            P.close()

        return figs


    def export_html_diags(self, htmlfile, **kwargs):
        """Make and export diagnostics as an html file

        Parameters
        ----------
        htmlfile: string
            Output html file name
        **kwargs
            All other params are passed to :meth:`plot_diags`

        Return
        ------
        string
            Html file name
        """
        # Get figures
        figs = self.plot_diags(**kwargs)
        figs = {'Ensemble diagnostics': figs}

        # Figure paths
        figs = dicttree_relpath(figs, os.path.dirname(htmlfile))

        # Render with template
        checkdir(htmlfile)
        render_and_export_html_template('dict2tree.html', htmlfile,
            title='Ensemble diagnostics', content=figs)
        self.created(htmlfile)
        return htmlfile

    def project_on_obs(self, obsmanager):
        """Interpolate the variables onto the observation locations or the
        :class:`~sonat.obs.ObsManager`
        """
        out = []
        for obs in obsmanager:
            vmod = [self.get_variable(vname, obs.depths)
                    for vname in obs.varnames]
            out.append(obs.project_model(vmod))
        return out


    def assert_compatible_with_obs(self, obsmanager, sync_norms=True,
            sync_missing_values=True):
        """Assert that the :class:`Ensemble` current instance is compatible
        with a :class:`~sonat.obs.ObsManager` instance

        It checks that observed variable are provided by the ensemble.
        It optinonally synchonise norms between model and observations.
        """

        # Check varnames
        for varname in obsmanager.suffixed_varnames:
            if varname not in self.varnames:
                raise SONATError(('Observed variable name "{}" not found '
                    ' in ensemble').format(varname))

        # Sync norms
        if sync_norms:
            self.sync_norms(force=False)
            dnorms = self.get_named_norms()
            obsmanager.set_named_norms(dnorms)

        # Missing values
        if sync_missing_values and self.missing_value != obsmanager.missing_value:
            obsmanager.missing_value = self.missing_value

        # Check projection
        oens = self.project_on_obs(obsmanager)
        poens = obsmanager.restack(oens)
        mask = N.isclose(poens, obsmanager.missing_value)
        if mask.ndim > 1:
            mask = mask.any(axis=1)
        valid = ~mask
        if not valid.any():
            raise SONATError("Ensemble does project at all onto observations")
        if not valid.all():
            self.warning("Ensemble does not fully projects onto observations")

        return dict(ens_on_obs=oens, packed_ens_on_obs=poens,
            packed_ens_on_obs_valid=valid)



