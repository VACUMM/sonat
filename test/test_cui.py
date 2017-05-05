"""Test script for module :mod:`sonat.ens`"""

from vcmq import netcdf4

from util import (CFGFILE)


from sonat.config import load_cfg
from sonat.cui import (
    ens_gen_pseudo_from_cfg,
    ens_plot_diags_from_cfg,
    obs_plot_from_cfg,
    arm_analysis_from_cfg,
    arm_sa_from_cfg,
    )

netcdf4()

CACHE = {}

def get_cfg():
    if 'cfg' in CACHE:
        return CACHE['cfg']
    CACHE['cfg'] = load_cfg(CFGFILE)
    return CACHE['cfg']


def test_cui_ens_gen_pseudo_from_cfg():

    # Load config
    cfg = get_cfg()

    # Run
    ens_gen_pseudo_from_cfg(cfg)

def test_cui_ens_plot_diags_from_cfg():

    # Load config
    cfg = get_cfg()

    # Run
    ens_plot_diags_from_cfg(cfg)

def test_cui_obs_plot_from_cfg():

    # Load config
    cfg = get_cfg()

    # Run
    obs_plot_from_cfg(cfg)

def test_cui_arm_analysis_from_cfg():

    # Load config
    cfg = get_cfg()

    # Run
    arm_analysis_from_cfg(cfg)

def test_cui_arm_sa_from_cfg():

    # Load config
    cfg = get_cfg()

    # Run
    arm_sa_from_cfg(cfg)

def _test_cui_test_from_cfg():
    from sonat import test_from_cfg

    # Load config
    cfg = get_cfg()

    # Run
    test_from_cfg(cfg)

if __name__=='__main__':

    test_cui_ens_gen_pseudo_from_cfg()
#    test_cui_ens_plot_diags_from_cfg()
#    test_cui_obs_plot_from_cfg()
#    test_cui_arm_analysis_from_cfg()
#    test_cui_arm_sa_from_cfg()
#    _test_cui_test_from_cfg()



