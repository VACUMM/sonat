"""The ARM interface"""

from .__init__ import SONATError
from .misc import _Base_
from .ens import Ensemble
from .obs import ObsManager, NcbobsPlatform
from ._fcore import f_arm

class ARM(_Base_):
    """Interface to the ARM assessment algorithm"""
    def __init__(self, ens, obsmanager, logger=None, **kwargs):

        # Init logger
        _Base_.__init__(self, logger=logger, **kwargs)

        # Check arguments
        if not isinstance(ens, Ensemble):
            raise SONATError('ens must be a Ensemble instance')
        msg = ('obsmanager must a either a ObsManager'
                    ' instance or a list of NcbobsPlatform instances')
        if isinstance(obsmanager, list):
            if not all([isinstance(obs, NcbobsPlatform) for obs in obsmanager]):
                raise SONATError(msg)
            obsmanager = ObsManager(obsmanager)
        elif not isinstance(obsmanager, ObsManager):
            raise SONATError(msg)
        self.ens = ens
        self.obsmanager = obsmanager

    def interp_ens_on_obs(self):
        """Interpolate the variables of an :class:`~sonat.ens.Ensemble` onto
        observation locations
        """
        out = []
        for obs in obsmanager:
            vmod = [ens.get_variable(vname, obs) for vname in obs.varnames]
            out.append(obs.interp_model(vmod))
        return out



    def get_matrices(self):
        """Get a dict of ensemble and observation matrices for ARM"""

        # Ensemble states on model grid
        Af = self.ens.stacked_data

        # Ensemble states at observation locations
        oens = self.interp_ens_on_obs()
        Yf = self.obsmanager.restack(oens)

        # Observation errors
        err = self.obsmanager.stacked_data
        R = N.asfortranarray(N.trace(err))

        return dict(Af=Af, Yf=Yf, R=R)


    def analyse(self):
        """Perform the ARM analysis"""

        # Get matrices
        mats = self.get_matrices()
        Af = mats['Af']
        Yf = mats['Yf']
        R = mats['R']

        # Call ARM
        ndof = min(dsamples.shape)
        spect, U, rho_mu, status = f_arm(ndof, Af, Yf, R)

