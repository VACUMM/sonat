"""The ARM interface"""
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
import inspect
import math
from collections import OrderedDict

import numpy as N
import MV2, cdms2

from vcmq import (curve, create_axis, dict_check_defaults, kwfilter)

from .__init__ import SONATError, sonat_warn
from .misc import _Base_, recursive_transform_att, vminmax
from .ens import Ensemble
from .obs import ObsManager, NcObsPlatform
from ._fcore import f_arm
from .pack import Simple1DPacker, default_missing_value

class ARM(_Base_):
    """Interface to the ARM assessment algorithm

    Example
    -------

    >>> arm = ARM(ens, obsmanager)
    >>> arm = ARM(ens, [ncobs1, ncobs2])

    >>> arm.analyse(ndof=5) # optional because on demand

    >>> print arm.nobs, arm.nstates, arm.nens, self.ndof
    >>> print arm.Af.shape # ensemble states
    >>> print arm.Yf.shape # ensemble states on obs
    >>> print arm.R.shape  # obs error

    >>> print arm.spect.shape # ARM spectrum
    >>> print arm.arm[0][0].shape # array modes for first platform and first var
    >>> print arm.rep[0] # representers for first var

    >>> print arm.raw_spect
    >>> print arm.raw_arm
    >>> print arm.raw_rep

    """
    def __init__(self, ens, obsmanager, logger=None, checkvars=True,
            syncnorms=True, missing_value=default_missing_value,
            bathy=None, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Check arguments
        if not isinstance(ens, Ensemble):
            raise SONATError('ens must be a Ensemble instance')
        msg = ('obsmanager must a either a ObsManager'
                    ' instance or a list of NcObsPlatform instances')
        if isinstance(obsmanager, list):
            if not all([isinstance(obs, NcObsPlatform) for obs in obsmanager]):
                raise SONATError(msg)
            obsmanager = ObsManager(obsmanager)
        elif isinstance(obsmanager, NcObsPlatform):
            obsmanager = ObsManager(obsmanager)
        elif not isinstance(obsmanager, ObsManager):
            raise SONATError(msg)
        self.ens = ens
        self.obsmanager = obsmanager

        # Sync missing values
        self._missing_value = missing_value
        self.set_missing_value(missing_value)

        # Check variables
        if checkvars:
            self.ens.check_variables()
            self.obsmanager.check_variables()

        # Check compatibility (names, norms, projection)
        diags = self.ens.assert_compatible_with_obs(self.obsmanager,
            syncnorms=syncnorms)
        self._ens_on_obs = diags['ens_on_obs']
        self._packed_ens_on_obs = diags['packed_ens_on_obs']
        self._packed_ens_on_obs_valid = diags['packed_ens_on_obs_valid']
        self._final_packer = Simple1DPacker(self._packed_ens_on_obs,
            self.missing_value, valid=diags['packed_ens_on_obs_valid'])

        # Inits
        self._inputs = {} # indirect input matrices
        self._results = {} # results

        # Bathy
        if bathy is not None:
            self.set_bathy(bathy)

    @property
    def nstate(self):
        return self.ens.stacked_data.shape[0]

    @property
    def nens(self):
        return self.ens.stacked_data.shape[1]

    @property
    def nobs(self):
        return self._final_packer.packed_data.shape[0]

    def set_missing_value(self, missing_value):
        if missing_value != self._missing_value:
            self.ens.missing_value = missing_value
            self.obsmanager.missing_value = missing_value

    def get_missing_value(self):
        return self._missing_value

    missing_value = property(fget=get_missing_value, fset=set_missing_value)
    fill_value = missing_value

    def set_bathy(self, bathy2d):
        """Interpolate a bathymetry and save it"""
        self.ens.set_bathy(bathy2d)
        self.obsmanager.set_bathy(bathy2d)

    def project_ens_on_obs(self):
        """Interpolate the variables of an :class:`~sonat.ens.Ensemble` onto
        observation locations
        """
        return self.ens.project_on_obs(self.obsmanager)

    @property
    def Af(self):
        """Ensemble states on model grid"""
        return self.ens.stacked_data

    @property
    def Yf(self):
        """Ensemble states at observation locations"""
        return self._final_packer.packed_data

    @property
    def R(self):
        """Observation errors"""
        if 'R' not in self._inputs:
            err = self._final_packer.pack(self.obsmanager.stacked_data)
            self._inputs['R'] = N.asfortranarray(N.diag(err))
        return self._inputs['R']

    @property
    def ndof(self):
        if 'raw_spect' in self._results:
            return len(self._results['raw_spect'])
        return min(self.Yf.shape)

    @property
    def mode_axis(self):
        if not hasattr(self, '_mode_axis') or len(self._mode_axis)!=self.ndof:
            self._mode_axis = create_axis((1, self.ndof+1), id='mode',
                                          long_name='Array mode id')
        return self._mode_axis

    def analyse(self, ndof=None):
        """Perform the ARM analysis

        Parameters
        ----------
        ndof: int, None
            Number of mode to retain. It defaults to the smallest dimension
            of the observation matrix.
        getraw: bool
            Also return raw array modes and array mode representers from ARM

        Return
        ------
        dict
            Keys are:

            - spect: spectrum
            - arm: array modes
            - rep: array mode representers
            - raw_arm: raw version of arm (optional)
            - raw_rep: raw version of rep (optional)

        """

        # Number of retained modes
        if ndof is None:
            ndof = self.ndof
        else:
            ndof = min(ndof, + self.ndof)

        # Check cache / truncate to ndof
        if 'spect' in self._results and ndof >= len(self._results[spect]):
            self._clean_analysis_(raw=False)
            self._results['spect'] = self._results['spect'][:ndof]
            self._results['arm'] = self._results['arm'][:, :ndof]
            self._results['rep'] = self._results['rep'][:, :ndof]
        else:
            self._clean_analysis_(raw=True)

        # Call ARM
        spect, arm, rep, status = f_arm(ndof, self.Af, self.Yf, self.R)

        # Check status
        if status:
            raise SONATError('Fortran ARM routine failed with exit status {}'
                .format(status))

        # Store raw results
        self._results = dict(
            raw_spect = spect,
            raw_arm = arm,
            raw_rep = rep,
        )

    def _clean_analysis_(self, raw):
        """Remove results from cache"""
        keys = ['spect', 'arm', 'rep']
        if raw:
            keys.extend(['raw_spect', 'raw_arm', 'raw_rep'])
        for key in keys:
            key = '_' + key
            if key in self._results:
                del self._results[key]

    def clean(self):
        self._clean_analysis_(True)

    def _check_analysis_(self):
        if 'raw_spect' not in self._results:
            self.analyse()

    @property
    def raw_spect(self):
        """Raw spectrum result"""
        self._check_analysis_()
        return self._results['raw_spect']

    @property
    def raw_arm(self):
        """Raw array mode results"""
        self._check_analysis_()
        return self._results['raw_arm']

    @property
    def raw_rep(self):
        """Raw representer results"""
        self._check_analysis_()
        return self._results['raw_rep']

    @property
    def spect(self):
        """ARM spectrum"""
        # Check cache
        if 'spect' in self._results:
            return self._results['spect']

        # Make sure that the analsys is done
        self._check_analysis_()

        # Unstack/pack/format
        self._results['spect'] = MV2.array(self.raw_spect,
            id='arm_spectrum', attributes=dict(long_name='ARM spectrum'),
            axes=[self.mode_axis])
        return self._results['spect']

    @property
    def arm(self):
        """Array modes"""
        # Check cache
        if 'arm' in self._results:
            return self._results['arm']

        # Make sure that the analsys is done
        self._check_analysis_()

        # Unstack/pack/format
        uarm = self._final_packer.unpack(self.raw_arm)
        self._results['arm'] = self.obsmanager.unstack(uarm, rescale=False,
            format=1, id='arm_{id}',
#            firstdims=[self.mode_axis],
            )
        def set_long_name(long_name):
            return 'Array modes of '+long_name
        recursive_transform_att(self._results['arm'], 'long_name', set_long_name)
        return self._results['arm']

    @property
    def rep(self):
        """Representers"""
        # Check cache
        if 'rep' in self._results:
            return self._results['rep']

        # Make sure that the analsys is done
        self._check_analysis_()

        # Unstack/pack/format
        self._results['rep'] = self.ens.unstack(self.raw_rep, rescale=False, format=1,
            id='rep_{id}',
#            firstdims=[self.mode_axis],
            )
        def set_long_name(long_name):
            return 'Array mode representers of '+long_name
        recursive_transform_att(self._results['rep'], 'long_name', set_long_name)
        return self._results['rep']

    def get_score(self, type='nev'):
        """Get the score for the current analysis

        Parameters
        ----------
        type: string
            Score type
       """
        return get_arm_score_function(type)(self.spect, self.arm, self.rep)

    def plot_spect(self, figfile='arm.spect.png', savefig=True, close=True,
                   title='ARM Sepctrum', hline=1., score='nev',
                   shade=True, **kwargs):
        """Plot the ARM spectrum"""
        kwhline = kwfilter(kwargs, 'hline')
        kwscore = kwfilter(kwargs, 'score')
        kwshade = kwfilter(kwargs, 'shade')

        # Main plot
        dict_check_defaults(kwargs, xmin=.5, xmax=self.ndof+.5, parg='o',
                            log=True, show=False, ymin=.01, fig='new')
        p = curve(self.spect, title=title,
              **kwargs)
        zorder = p['curve'][0].zorder

        # Horizontal line at 1
        if hline:

            # Shading
            if shade:
                dict_check_defaults(kwshade, linewidth=0, color='.8',
                                    zorder=zorder-0.01)
                axis = p.axes.axis()
                p.axes.axhspan(axis[2], hline, **kwshade)
                p.axes.axis(axis)

            # Line
            dict_check_defaults(kwhline, linewidth=.5, color='k',
                                zorder=zorder+0.01)
            p.axes.axhline(1., **kwhline)

        # Score
        if score:
            score = self.get_score(score)
            dict_check_defaults(kwscore, family='monospace', ha='right', va='top',
                                bbox=dict(linewidth=0, facecolor='w', alpha=.2))
            p.text(.95, .95, 'Score: {:.01f}'.format(score), transform='axes',
                   **kwscore)

        # Save
        if savefig:
            p.savefig(figfile)
            self.created(figfile)
            return {'Spectrum':figfile}


    def plot_arm(self, varnames=None,
                 figpat='arm.arm.mode{mode}_{var_name}_{slice_type}_{slice_loc}.png',
                 **kwargs):
        """Plot the array modes"""
        self.verbose('Plotting array modes')

        # Select data
        variables = self.obsmanager.select_variables(varnames, source=self.arm,
                                                     prefix_to_rm='arm_')

        # Plot
        modemax = self.format_mode()
        figs = OrderedDict()
        subfigs = figs['Array modes'] = OrderedDict()
        dict_check_defaults(kwargs, sync_vminmax='symetric', cmap='cmocean_balance')
        for imode in range(self.ndof):
            mode = self.format_mode(imode)
            self.debug(' Plotting array mode {}/{}'.format(mode, modemax))

            # Extract single mode
            kwargs['subst'] = {'mode':mode, 'modemax':modemax}
            mvars = self.extract_mode(variables, imode)
            remove_arm = lambda id: id[4:]
            recursive_transform_att(mvars, 'id', remove_arm)

            # Make plots
            mfigs = self.obsmanager.plot(mvars, input_mode='arrays',
                                         title='Array mode {mode}/{modemax} - {var_name} - {slice_loc}',
                                         figpat=figpat, **kwargs)

            # Register figures
            if isinstance(mfigs, dict):
                subfigs['Mode '+mode] = mfigs

        return figs

    def plot_rep(self, varnames=None, imodes=None,
                 add_obs=True,
                 titlepat='Representer {mode}/{mode_max} - {var_name} - {slice_loc}',
                 figpat='arm.rep.mode{mode}_{var_name}_{slice_type}_{slice_loc}.png',
                 sync_vminmax=True,
                 **kwargs):
        """Plot the representers of array modes

        Parametes
        ---------
        varnames: strings
            List variable names to select
        add_obs: bool
            Add the location of all observations?
        figpat: string
            Output figure file pattern
        sync_vminmax: bool
            Sync min and max for all plots
        """
        self.verbose('Plotting representers of array modes')

        # Select data
        variables = self.ens.select_variables(varnames, source=self.rep,
                                              prefix_to_rm='rep_')

        # Plot
        mode_max = self.format_mode()
        figs = OrderedDict()
        subfigs = figs['Representer of array modes'] = OrderedDict()
        dict_check_defaults(kwargs, cmap='cmocean_balance',
                            levels_mode='symetric')
        if sync_vminmax:
             dict_check_defaults(kwargs, **vminmax(self.rep, asdict=True))
        if imodes is None:
            imodes = range(self.ndof)
        elif N.isscalar(imodes):
            imodes = [imodes]
        for imode in imodes:

            mode = self.format_mode(imode)
            self.debug(' Plotting representer of mode {}/{}'.format(mode, mode_max))

            # Extract single mode
            mvars = self.extract_mode(variables, imode)

            # Sliced plot
            kwargs['subst'] = {'mode':mode, 'mode_max':mode_max}
            if add_obs:
                kwargs['obs'] = (self.obsmanager
                       if isinstance(add_obs, bool) else self.obsmanager[add_obs])
            ff = self.ens.plot_fields(mvars, figpat=figpat, titlepat=titlepat,
                                      fig='new', **kwargs)

            # Register figures
            if isinstance(ff, dict):
                subfigs['Mode '+mode] = ff

        return figs




    def format_mode(self, imode=None, fmt='{imode:0{ndigit}d}'):
        if imode is None:
            imode = self.ndof-1
        ndigit = int(math.log10(self.ndof)) + 1
        return fmt.format(**locals())


    def extract_mode(self, data, imode):
        """Extract a single mode from recursive lists of arrays

        The mode axis is supposed to be the first one.
        """
        if not isinstance(data, list):
            dat = data[imode]
            return dat
        return [self.extract_mode(dat, imode) for dat in data]


#: Prefix of all ARM score functions
ARM_SCORE_FUNCTION_PREFIX = 'arm_score_'

def arm_score_nev(spect, arm, rep):
    """ARM score as the number of eigen values greater than one"""
    return (spect>=1).sum()

def arm_score_fnev(spect, arm, rep):
    """ARM score as the fractional number of eigen values greater than one"""
    # Integer value
    nev = arm_score_nev(spect, arm, rep)

    # Linear interpolation
    if nev==0 or nev==len(spect):
        return float(nev)
    return float(nev) + spect[nev-1] / (spect[nev-1] - spect[nev])


def arm_score_relvar(spect, arm, rep):
    """ARM score as the fractional number of eigen values greater than one"""
    # Float value
    fnev = arm_score_fnev(spect, arm, rep)
    if fnev==0:
        return 0.
    nev = int(fnev)

    # One-based spectum
    spect = spect - 1

    # Base
    var = spect[:nev].sum()
    var -= spect[0] * .5
    var -= spect[-1] * .5
    if int(fnev)==len(spect):
        return var

    # Fraction
    var += fnev//1 * spect[-1] * .5

#: Registered ARM score functions
ARM_SCORE_FUNCTIONS = {}

def register_arm_score_function(func, warn=True, replace=True):
    """Register a new score function"""

    # Get the generic name
    fname = func.__name__.lower()
    prefix = 'arm_score_'
    if fname.startswith(prefix):
        fname = fname[len(prefix):]

    # Inspect
    argspec = inspect.getargspec(func)
    nargs = len(argspec[0])
    if argspec[3] is not None:
        nargs -= len(argspec[3])
    if nargs != 3:
        raise SONATError("You score function must exactly 3 positional "
            "arguments,  not {}: spect, arm, rep".format(nargs))

    # Save it
    if fname in ARM_SCORE_FUNCTIONS:
        if warn:
            msg = 'ARM score function already registered: {}'.format(fname)
            if replace:
                msg = msg + '. Replacing it...'
            sonat_warn(msg)
        if not replace:
            return
    ARM_SCORE_FUNCTIONS[fname] = func

    return fname, func

def get_arm_score_function(fname):
    """Get a score function from its name"""

    fname = fname.lower()
    if fname not in ARM_SCORE_FUNCTIONS:
        raise SONATError('Unkown score function "{}". Valid names: {}'.format(
            fname, ', '.join(ARM_SCORE_FUNCTIONS.keys())))

    return ARM_SCORE_FUNCTIONS[fname]

# Register default score functions
register_arm_score_function(arm_score_nev)
register_arm_score_function(arm_score_fnev)
register_arm_score_function(arm_score_relvar)
