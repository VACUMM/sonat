"""The ARM interface"""
import numpy as N
import MV2

from .__init__ import SONATError
from .misc import _Base_
from .ens import Ensemble
from .obs import ObsManager, NcObsPlatform
from ._fcore import f_arm

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
            syncnorms=True, **kwargs):

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

        # Check variables
        if checkvars:
            self.ens.check_variables()
            self.obsmanager.check_variables()

        # Make sure that ens norms are synced with obs norms
        if syncnorms:
            self.ens.sync_norms(force=False)
            dnorms = self.ens.get_named_norms()
            self.obsmanager.set_named_norms(dnorms)

        # Dims
        self.nstate, self.nens = self.ens.stacked_data.shape
        self.nobs = self.obsmanager.stacked_data.shape[0]

        # Inits
        self._inputs = {} # indirect input matrices
        self._results = {} # results

    def project_ens_on_obs(self):
        """Interpolate the variables of an :class:`~sonat.ens.Ensemble` onto
        observation locations
        """
        out = []
        for obs in self.obsmanager:
            vmod = [self.ens.get_variable(vname, obs) for vname in obs.varnames]
            out.append(obs.project_model(vmod))
        return out

    @property
    def Af(self):
        """Ensemble states on model grid"""
        return self.ens.stacked_data

    @property
    def Yf(self):
        """Ensemble states at observation locations"""
        if 'Yf' not in self._inputs:
            oens = self.project_ens_on_obs()
            self._inputs['Yf'] = self.obsmanager.restack(oens)
        return self._inputs['Yf']

    @property
    def R(self):
        """Observation errors"""
        if 'R' not in self._inputs:
            err = self.obsmanager.stacked_data
            self._inputs['R'] = N.asfortranarray(N.diag(err))
        return self._inputs['R']

    @property
    def ndof(self):
        if 'spect' in self._results:
            return len(self._results[spect])
        return min(self.Yf.shape)


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
            self.analysis()

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
        self._results['spect'] = MV2.array(spect, id='arm_spectrum',
            attributes=dict(long_name='ARM spectrum'))
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
        self._results['arm'] = self.obsmanager.unstack(arm, rescale=False,
            format=1, id='arm_{id}')
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
        self._results['rep'] = self.ens.unstack(rep, rescale=False, format=1,
            id='arm_rep_{id}')
        return self._results['rep']

    def get_score(self, type='nev'):
        """Get the score for the current analysis

        Parameters
        ----------
        type: string
            Score type
       """
        return get_arm_score_function(type)(self.spect, self.arm, self.rep)



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


ARM_SCORE_FUNCTIONS = {}

def register_arm_score_function(func):
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
        sonat_warn('ARM score function "{}" is already registered. Replacing...'.format(
            fname))
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
