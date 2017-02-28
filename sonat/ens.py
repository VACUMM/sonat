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
from vcmq import (cdms2, MV2, DS, ncget_time, lindates, ArgList, format_var,
    MV2_concatenate, create_axis, N, kwfilter, check_case, match_known_var,
    CaseChecker, curve, bar, dict_check_defaults, dict_merge, checkdir,
    dicttree_get)

from .__init__ import get_logger, sonat_warn, SONATError
from .misc import (list_files_from_pattern, ncfiles_time_indices, asma, NcReader,
    validate_varnames, _NamedVariables_, check_variables, dicttree_relpath,
    interpret_level)
from .stack import Stacker
from ._fcore import f_eofcovar, f_sampleens
from .plot import (plot_gridded_var, DEFAULT_PLOT_KWARGS)
from .render import register_html_template, render_and_export_html_template

#: Default plotting keyword per ensemble diagnostic
DEFAULT_PLOT_KWARGS_PER_ENS_DIAG = {
    'explained_variance':dict(
        xmin=.5, ymin=0,
        title='Explained variance', width=.5,
        xlabel="%(long_name)s",
        bottom=.2,
    ),
    'local_explained_variance':{
        'cmap':'speed',
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
    depths: string, None, list of floats, array, tuple of them, dict
        Here are some possible values:

        - "surf" or "bottom": self explenatory
        - None or "3d": No slice, so get all levels with no interpolation.
        - A list or array of negative depths: get all levels and
          interpolate at these depths.

        Variables sliced with "surf" and "bottom" are returned with
        an id suffixed with "_surf" or "_bottom".
        You can speficy different slicing  using a tuple
        of depths specifications.
        You can specilise slicings a variable using a dictionary with
        the key as the variable name.



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
    syncnorms: bool, optional
        Make sur that all variable whose id has the same prefix (string
        before "_") have the same norm
    """

    def __init__(self, data, ev=None, variance=None, logger=None, norms=None,
            means=None, syncnorms=True, checkvars=False, **kwargs):

        # Check input variables
        datas = ArgList(data).get()
        if checkvars:
            check_variables(datas)
        for var in datas: # rename first axis
            var.getAxis(0).id = 'member'

        # Init base
        kwargs['nordim'] = False
        Stacker.__init__(self, data, logger=logger,
                         norms=None if isinstance(norms, dict) else norms,
                         means=means, **kwargs)
        self.variables = self.inputs

        # Named norms
        if norms and isinstance(norms, dict):
            self.set_named_norms(norms)

        # Synchronise norms
        self._norm_synced = False
        if syncnorms:
            self.sync_norms()

        # Add more variables
        self.ev = ev
        self.variance = variance

        # Caches
        self._diags = {}
        self._meananoms = {}


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
        diag_dict[diag_name] = self.unmap(dds)


    def get_diags(self, mean=True, variance=True, kurtosis=True, skew=True,
            skewtest=True, kurtosistest=True, normaltest=True):
        """Get ensemble diagnostics as a dict

        Example
        -------
        >>> diags = ens.get_diags(mean=False, skew=True)
        >>> assert diags['skew'].shape == diags['skewtest'].shape

        Sea also
        --------
        :func:`scipy.stats.kurtosis` :func:`scipy.stats.kurtosistest`
        :func:`scipy.stats.skew` :func:`scipy.stats.skewtest`
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
                diags['variance'] = self._diags['variance'] = self.unmap(vars)
                diags['vars'] = vars

            # Explained variance
            if self.ev is not None and hasattr(self.ev, 'total_variance'):
                if 'explained_variance' in self._diags:
                    diags['explained_variance'] = self._diags['explained_variance']
                else:
                    diags['explained_variance'] = self._diags['explained_variance'] = \
                        100 * (self.ev**2) / self.ev.total_variance
                    diags['explained_variance'].id = 'explained_variance'
                    diags['explained_variance'].units = '%'

            # Local explained variance
            if self.variance is not None:
                if 'local_explained_variance' in self._diags:
                    diags['local_explained_variance'] = self._diags['local_explained_variance']
                else:
                    evars = []
                    for lvar, var in zip(diags['vars'], self.remap(self.variance)):
                        if cdms2.isVariable(lvar):
                            evar = lvar.clone()
                            evar.id = var.id + '_variance'
                            evar.units = '%'
                        else:
                            evar = lvar.copy()
                        evar[:] /= var
                        evar[:] *= 100
                        evars.append(evar)
                    diags['local_explained_variance'] = self._diags['local_explained_variance'] = \
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
            titlepat = '%(long_name)s - {diag_longname} - {slice_loc}',
            depths=None, points=None,
            zonal_sections=None, merid_sections=None,
            figpat_slice='sonat.ens.{diag_name}_{var_name}_{slice_type}_{slice_loc}.png',
            figpat_generic='sonat.ens.{diag_name}.png',
            fmtlonlat='{:.2f}', fmtdep='{:04.0f}m',
            figfir=None, show=False,
            htmlfile=None, props=None):
        """Create figures for diagnostics

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

        # Diags
        diags = self.get_diags(mean=mean, variance=variance, kurtosis=kurtosis, skew=skew,
            skewtest=skewtest, kurtosistest=kurtosistest, normaltest=normaltest)

        # Inits
        self.verbose('Plotting ensemble diagnostics')
        figs = OrderedDict()
        if depths and not isinstance(depths, list):
            depths = [depths]
        if props is None:
            props = {}
        kwmap = props.get('map', {})
        kwcurve = props.get('curve', {})
        kwprops = dict_merge(props, DEFAULT_PLOT_KWARGS_PER_ENS_DIAG)
        depths = ArgList(depths).get()
        depths3d = [dd for dd in depths if dd not in ['surf', 'bottom']]
        depths3d.sort(reverse=True)

        # Loop on diags
        figs = OrderedDict()
        for diag_name, diag_var in diags.items():

            diag_longname = diag_name.title().replace('_', ' ')
            self.debug('Plotting ens diag: '+diag_longname)

            # Explained variance
            if diag_name=='explained_variance':

                figfile = figpat_generic.format(**locals())
                checkdir(figfile)
                kw = kwprops.get(diag_name, {})
                dict_check_defaults(kw, xmax=len(diag_var)+.5, savefig=figfile,
                    **DEFAULT_PLOT_KWARGS)
                bar(diag_var, **kw)
                self.created(figfile)
                figs[diag_longname] = figfile

            # Loop on variables
            else:

                for diag_var in ArgList(diag_var).get():
                    order = diag_var.getOrder()

                    var_name, vdepth, vdiag_name = split_varname(diag_var)
                    varname = var_name
                    id = diag_var.id
                    self.debug(' Variable: '+var_name)


                    # Maps
                    if depths is not False:
                        toplot = []
                        slice_type = 'map'

                        # Get what to plot
                        if ((vdepth=='surf' and 'surf' in depths) or
                                (vdepth=='bottom' and 'bottom' in depths)): # 2D

                            toplot.append(dict(loc=vdepth, var=diag_var))

                        elif diag_var.getLevel() is not None and depths3d: # 3D

                            # Loop on 3D depths specs
                            for d3d in depths3d:

                                if d3d in [None, '3d']: # simple slices

                                    levels = diag_var.getLevel()[::-1]
                                    for i, level in enumerate(levels):
                                        toplot.append(dict(
                                            loc=fmtdep.format(abs(level)),
                                            var=diag_var[len(levels)-i-1]))

                                else: # interpolations

                                    levels = N.sort(N.atleast_1d(d3d))
                                    for i, level in enumerate(levels):
                                        toplot.append(dict(
                                            loc=fmtdep.format(abs(level)),
                                            var=diag_var,
                                            kw=dict(depth=level)))

                        # Init dict
                        if toplot:
                            dd = figs.setdefault(
                                diag_longname, OrderedDict()).setdefault(
                                    var_name.upper(),  OrderedDict()).setdefault(
                                        slice_type.title(), OrderedDict())

                        # Plot them
                        for spec in toplot:
                            slice_loc = spec['loc']
                            var = spec['var']
                            self.debug('  Location: '+slice_loc)
                            title = titlepat.format(**locals())
                            figfile = figpat_slice.format(**locals())
                            checkdir(figfile)
                            kw = kwprops.get(diag_name, {})
                            if 'kw' in spec:
                                kw.update(**spec['kw'])
                            kw.update(title=title, savefig=figfile, **kwmap)
                            dict_check_defaults(kw, **DEFAULT_PLOT_KWARGS_PER_ENS_DIAG)
                            plot_gridded_var(diag_var, **kw)
                            self.created(figfile)
                            dd[slice_loc] = figfile
                        del toplot

        # Export to html
        if htmlfile:

            # Figure paths
            figs = dicttree_relpath(figs, os.path.dirname(htmlfile))

            # Render with template
            checkdir(htmlfile)
            render_and_export_html_template('dict2tree.html', htmlfile,
                title="Ensemble diagnostics", content=figs)
            self.created(htmlfile)

        return figs

    def export_diags(self, htmlfile, **kwargs):
        """Make and export diagnostics as an html file

        Parameters
        ----------
        htmlfile: string
            Output html file name
        **kwargs
            All other params are passed to :meth:`plot_diags`

        """
        kwargs['htmlfile'] = htmlfile
        figs = self.plot_diags(**kwargs)
        return htmlfile, figs

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


    def assert_compatible_with_obs(self, obsmanager, syncnorms=True,
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
        if syncnorms:
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

