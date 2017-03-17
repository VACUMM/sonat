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
import re
from collections import OrderedDict

import numpy as N
import MV2, cdms2

from vcmq import (curve, create_axis, dict_check_defaults, kwfilter, m2deg)

from .__init__ import SONATError, sonat_warn
from .misc import _Base_, recursive_transform_att, vminmax, xycompress
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
            sync_norms=True, missing_value=default_missing_value,
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
        self.obsmanager = self.obs = obsmanager

        # Sync missing values
        self._missing_value = missing_value
        self.set_missing_value(missing_value)

        # Check variables
        if checkvars:
            self.ens.check_variables()
            self.obsmanager.check_variables()

        # Check compatibility (names, norms, projection) and project
        self.project_ens_on_obs(sync_norms=sync_norms)

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

    def project_ens_on_obs(self, sync_norms=True, sync_missing_values=True):
        """Interpolate the variables of an :class:`~sonat.ens.Ensemble` onto
        observation locations
        """
        diags  = self.ens.assert_compatible_with_obs(self.obsmanager,
            sync_norms=sync_norms, sync_missing_values=sync_missing_values)
        self._ens_on_obs = diags['ens_on_obs']
        self._packed_ens_on_obs = diags['packed_ens_on_obs']
        self._packed_ens_on_obs_valid = diags['packed_ens_on_obs_valid']
        self._final_packer = Simple1DPacker(diags['packed_ens_on_obs'],
            self.missing_value, valid=diags['packed_ens_on_obs_valid'])
#        return self.ens.project_on_obs(self.obsmanager)

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
    def S(self):
        if 'S' not in self._inputs:
            err = self._final_packer.pack(self.obsmanager.stacked_data)
            Rinv = N.diag(1/err)
            self._inputs['S'] = N.dot(N.sqrt(Rinv / (self.ndof-1)), self.Yf)
        return self._inputs['S']

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

    def analyse(self, ndof=None, force=False):
        """Perform the ARM analysis

        Parameters
        ----------
        ndof: int, None
            Number of mode to retain. It defaults to the smallest dimension
            of the observation matrix.
        force: bool
            Force the analysis without checking cache

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
        if (not force  and 'raw_spect' in self._results and
            ndof >= len(self._results[raw_spect])):
            self.clean_analysis(raw=False)
            self._results['raw_spect'] = self._results['raw_spect'][:ndof]
            self._results['raw_arm'] = self._results['raw_arm'][:, :ndof]
            self._results['raw_rep'] = self._results['raw_rep'][:, :ndof]
            return
        else:
            self.clean_analysis(raw=True)

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

    def clean_analysis(self, raw):
        """Remove results from cache"""
        keys = ['spect', 'arm', 'rep']
        if raw:
            keys.extend(['raw_spect', 'raw_arm', 'raw_rep'])
        for key in keys:
            if key in self._results:
                del self._results[key]

    def clean_inputs(self):
        self._inputs = {}

    def clean(self):
        self.clean_analysis(True)
        self.clean_inputs()

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

#: Sensitivy analyser
class ARMSA(_Base_):

    def __init__(self, arm, logger=None, **kwargs):

        # Check arm
        if not isinstance(arm, ARM):
            raise SONATError('Input must be an ARM instance')

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        self.arm = arm
        self.obs = arm.obs
        self.ens = arm.ens


    def export(self, **kwargs):

        raise SONATError('The method must be overwritten')

RE_PERTDIR_MATCH = re.compile(r'[+\-][xy]$').match

class XYARMSA(ARMSA):


    def __init__(self, arm, **kwargs):

        ARMSA.__init__(self, arm, **kwargs)

        if not any([obs.mobile for obs in self.arm.obs]):
            self.warning("None of the platforms is mobile")

        # Backup stuff
        self.backups = {}
        # - platform coordinates
        self.backups['obs'] = {}
        for obs in self.arm.obs:
            if obs.mobile:
                self.backups['obs'][obs] = {}
                for att in 'lons', 'lats':
                    self.backups['obs'][obs][att] = getattr(obs, att).copy()
        # - ARM
        self.backups['arm'] = {}
        for att in ('_ens_on_obs', '_packed_ens_on_obs',
                    '_packed_ens_on_obs_valid', '_final_packer'):
            if hasattr(self.arm, att):
                self.backups['arm'][att] = getattr(self.arm, att)

    def restore_arm(self):
        for att, val in self.backups['arm'].items():
            setattr(self.arm, att, val)

    def restore_platform(self, obs):
        if not obs in self.backups['obs']:
            return
        for att, val in self.backups['obs'][obs].items():
            setattr(obs, att, val)

    def restore(self):
        self.restore_arm()
        for obs in self.obs:
            self.restore_platform(obs)

    def activate_platform_pert(self, obs, pert, direction, index):
        """Change the horizontal location of the platform at a single point
        and in a single direction.

        Parameters
        ----------
        pert: float
            Change in position in meridional degrees
        direction: string
            Possible values:

            - "[+|-][x|y]": negative or positive perturbation along x or y
            - "": back to original position

        index: int
            Index of a mobile observation.
            Its coordinates are accessed at ``obs.lons[index]`` and
            ``obs.lats[index]``.
        """
        if not obs.mobile:
            return

        if direction and pert:

            # Parse direction
            direction = str(direction).lower()
            if not RE_PERTDIR_MATCH(direction):
                raise SONATError("The direction of perturbation argument must"
                    " be of the form: '{+|-}{x|y}'")
            isx = direction[1]=='x'
            sign = 1 if direction[0]=='-' else -1

            # Reset
            self.deactivate_platform_pert(obs)

            # X and Y perturbation values
            if isx:
                lat = obs.lats.mean()
                pert = m2deg(pert, lat)
                cname = 'lons'
            else:
                cname = 'lats'
            coords = getattr(obs, cname) # array of coordinates

            # Change coordinates in place
            coords[index] += pert

        else: # back to orig

            self.deactivate_platform_pert()


    def deactivate_platform_pert(self, obs):
        """Get back to original locations"""
        self.restore_platform(obs)

    def get_platform_pert_indices_iter(self, obs):
        """Get an iterator on mobile indices

        ``obs.lons[index]`` and ``obs.lats[index]`` are the coordinates
        of a mobile point.
        """
        if not obs.mobile:
            return iter([])
        if not hasattr(obs, '_xy_pert_indices'):
            obs._xy_pert_indices = N.where(obs.mobility>=1)[0]
        return iter(obs._xy_pert_indices)

    def init_platform_pert_arrays(self, obs):
        """Init the arrays that will filled by the sensitivity analysis

        Return
        ------
        array
            Northward derivative
        array
            Southward derivative
        array
            Eastward derivative
        array
            Westward derivative
        """
        if not obs.mobile:
            sonat_warn('This platform cannot be moved')
            return

        mobile = obs.mobility>=1
        dsn = obs.mobility.clone().astype('d')
        dsn[:] = N.ma.masked
        dsn.id = 'dsn'
        dsn.long_name = 'Northward derivative of score'
        dss = dsn.clone()
        dss.id = 'dss'
        dss.long_name = 'Southward derivative of score'
        dsw = dsn.clone()
        dsw.id = 'dsw'
        dsw.long_name = 'Westward derivative of score'
        dse= dsn.clone()
        dse.id = 'dse'
        dse.long_name = 'Eastward derivative of score'

        return dse, dsw, dsn, dss


    def analyse(self, platforms=None, pert=0.001, score_type='fnev', direct=False):
        """Perform the senisitivty analysis

        Parameters
        ----------
        platforms: strings
            Selected platform names
        pert: float
            Location perturbation in meridional degrees
        score_type: string
            One of the registered score types
        direct: bool
            If True, a full ARM analysis is made after each perturbation,
            otherwise the spectrum is computed by projecting the new S matrix
            onto the original EOFs, which is faster.

        Return
        ------
        dict
            Keys are the mobile platforms and values are the four sensitivity
            arrays in each direction: east, west, north, south.
            Each array has the size of the number of platform mobile locations.
        """


        # Reference
        score0 = self.arm.get_score(score_type)
        ev0=self.arm.raw_spect
        r0=self.arm.R.copy()
        Yf0=self.arm.Yf.copy()

        # Loop on platforms
        results = {}
        score_func = get_arm_score_function(score_type)
        for obs in self.obs:

            # Nothing to move
            if not obs.mobile:
                continue

            # Selected platforms only
            if platforms is not None:
                if isinstance(platforms, str):
                    platforms = [platforms]
                if not obs.name in platforms:
                    continue

            # Inits
            dse, dsw, dsn, dss = self.init_platform_pert_arrays(obs)
            results[obs] = dse, dsw, dsn, dss

            # Loop on mobile points
            for idx in self.get_platform_pert_indices_iter(obs):

                # Check mobility
                if not obs.mobility[idx]:
                    continue

                # Loop on directions
                for ds, pdir in zip([dse, dsw, dsn, dss],
                                    ['+x', '-x', '+y', '-y']):

                    # Change location
                    self.activate_platform_pert(obs, pert, pdir, idx)

                    # Project ensemble on obs
                    self.arm.project_ens_on_obs(sync_norms=False,
                                                sync_missing_values=False)
                    #FIXME: intercept changes in obs shape

                    # Method
                    if direct:

                        # Reset ARM analysis
                        self.arm.clean()

                        # Analyse and get score
                        score1 = self.arm.get_score(score_type)
                        ev1=self.arm.raw_spect
                        r1=self.arm.R.copy()
                        Yf1=self.arm.Yf.copy()

                    else:

                        # Clean inputs only
                        self.arm.clean_inputs()

                        # Prevent change of size errors due to change of locations
                        try:

                            # Projection onto EOFs
                            pcs = N.dot(self.arm.S.T, self.arm.raw_arm)

                            # Indirect spectrum
                            spect = (pcs**2).sum(axis=0)

                            # Score
                            score1 = score_func(spect, self.arm.raw_arm,
                                                self.arm.raw_rep)

                        except:

                            score1 = N.ma.masked


                    # Fill array
                    print score1
                    ds[idx] = (score1 - score0) / pert
                    print ds[idx]

        self.restore()

        return results

    def plot(self, platforms=None, pert=0.01, score_type='fnev'):

        results = self.analysis(platforms=platforms, pert=pert,
                                score_type=score_type, method=method)



#: Registered ARM sensitivity analysers
ARM_SENSITIVITY_ANALYSERS = {}

def register_arm_sensitivity_analyser(cls, name=None, warn=True, replace=False):

    # Check
    if not issubclass(cls, ARMSA):
        raise SONATError('cls must be a subclass of ARMSA')

    # Name
    if name is None:
        name = cls.__class__.__name__.lower()
        if name.endswith('armsa'):
            name = name[:-5]

    # Register
    if name in ARM_SENSITIVITY_ANALYSERS:
        if warn:
            msg = 'ARM sensitivity analyser already registered: {}'.format(name)
            if replace:
                msg = msg + '. Replacing it...'
            sonat_warn(msg)
        if not replace:
            return
    ARM_SENSITIVITY_ANALYSERS[name] = cls
    cls.name = name
    return name, cls

register_arm_sensitivity_analyser(XYARMSA)

def get_arm_sensitivity_analyser(self, name, arm=None):
    """Get an ARM sensitivity analyser class or instance"""
    if name not in ARM_SENSITIVITY_ANALYSERS:
        raise SONATerror(('Invalid ARM sensitivity analyser: {}. '
                          'Please choose one of: {}').format(name,
                          ' '.join(ARM_SENSITIVITY_ANALYSERS.keys())))

    # Class
    cls = ARM_SENSITIVITY_ANALYSERS[name]
    if arm is None:
        return cls

    # Instance
    return cls(arm)

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
