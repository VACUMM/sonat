"""The ARM interface"""
import numpy as N
import MV2

from .__init__ import SONATError
from .misc import _Base_
from .ens import Ensemble
from .obs import ObsManager, NcObsPlatform
from ._fcore import f_arm

class ARM(_Base_):
    """Interface to the ARM assessment algorithm"""
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
        self.nobs = self.ens.stacked_data.shape[0]

    def project_ens_on_obs(self):
        """Interpolate the variables of an :class:`~sonat.ens.Ensemble` onto
        observation locations
        """
        out = []
        for obs in self.obsmanager:
            vmod = [self.ens.get_variable(vname, obs) for vname in obs.varnames]
            out.append(obs.project_model(vmod))
        return out



    def get_matrices(self):
        """Get a dict of ensemble and observation matrices for ARM"""

        # Ensemble states on model grid
        Af = self.ens.stacked_data

        # Ensemble states at observation locations
        oens = self.project_ens_on_obs()
        Yf = self.obsmanager.restack(oens)

        # Observation errors
        err = self.obsmanager.stacked_data
        R = N.asfortranarray(N.diag(err))

        return dict(Af=Af, Yf=Yf, R=R)


    def analyse(self, getraw=False):
        """Perform the ARM analysis

        Parameters
        ----------
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

        # Get matrices
        mats = self.get_matrices()
        Af = mats['Af']
        Yf = mats['Yf']
        R = mats['R']

        # Call ARM
        ndof = min(Yf.shape)
        spect, arm, rep, status = f_arm(ndof, Af, Yf, R)

        # Untack/pack/format to output
        out = {}
        # - spectrum
        out['spect'] = MV2.array(spect, id='arm_spectrum', attributes=dict(
            long_name='ARM spectrum'))
        # - array modes
        out['arm'] = self.obsmanager.unstack(arm, rescale=False, format=1,
            id='arm_{id}')
        if getraw:
            out['raw_arm'] = arm
        # - representers
        out['rep'] = self.ens.unstack(rep, rescale=False, format=1,
            id='arm_rep_{id}')
        if getraw:
            out['raw_rep'] = rep

        return out
