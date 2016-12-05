"""Test script for module :mod:`sonat.ens`"""

import os
import sys
import numpy as N
import cdms2
import cdtime
from vcmq import comptime, netcdf4

from util import (assert_allclose, LOGGER, NCFILE_OBS_SURF)

from sonat.obs import NcObsPlatform


def test_obs_ncobsplatform():

    # Load
    nop = NcObsPlatform(NCFILE_OBS_SURF, logger=LOGGER, lat=(45, 47.8))

    assert nop.lons.shape == (3, )

    return nop

if __name__=='__main__':
    nop = test_obs_ncobsplatform()
