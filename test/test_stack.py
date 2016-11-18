"""Test script for module :mod:`pyarm.stack`
"""
import os
import sys
import numpy as N
import MV2

assert_allclose = N.testing.assert_allclose

thisdir = os.path.dirname(__file__)
libdir = os.path.abspath(os.path.join('..', 'pyarm'))
if os.path.exists(os.path.join(libdir, '__init__.py')):
    sys.path.insert(0, os.path.dirname(libdir))
from pyarm import get_logger
from pyarm.stack import Stacker
from util import create_mv2_gridder_xyzt, create_mv2_scattered_xyzt

logger = get_logger(level="debug")

def test_stack_mv2_with_record():

    # Fake data
    # - first array
    data0 = create_mv2_gridder_xyzt(rotate=30)
    data0[:, :, 3:5, 2:4] = MV2.masked
    raxis = data0.getTime()
    del raxis.axis
    del raxis.units
    raxis.id = 'member'
    # - second array
    data1 = create_mv2_gridder_xyzt(rotate=0, nx=10, ny=9, nz=4)
    data1[:, :, 5:7, 6:7] = MV2.masked
    data1.setAxis(0, raxis)

    # Stack
    stacker = Stacker([data0, data1], nordim=False, logger=logger)

    # Unstack
    unstacked0, unstacked1 = stacker.unstack(stacker.stacked_data)

    # Restack
    restacked = stacker.restack([data0, data1])

    # Checks
    assert_allclose(stacker.stacked_data.shape,
        ((~data0[0].mask).sum()+(~data1[0].mask).sum(), data0.shape[0]))
    assert_allclose(unstacked0, data0)
    assert_allclose(unstacked1, data1)
    assert_allclose(restacked, stacker.stacked_data)

    return stacker

def test_pack_mv2_scattered_without_record_fixed_norm():

    # Fake data
    # - first array
    lons, lats, data = create_mv2_scattered_xyzt(nt=0)
    data[:, 3:5] = MV2.masked
    data[:] *= data
    norm = data.std()*2.

if __name__=='__main__':
    stacker = test_stack_mv2_with_record()
    print 'Done'
