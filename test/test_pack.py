"""Test script for module :mod:`pyarm.pack`
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
from pyarm.pack import Packer
from util import create_mv2_gridder_xyzt, create_mv2_scattered_xyzt

logger = get_logger(level="debug")

def test_pack_mv2_curved_with_record():

    # Fake data
    data = create_mv2_gridder_xyzt(rotate=30)
    data[:, :, 3:5, 2:4] = MV2.masked
    raxis = data.getTime()
    del raxis.axis
    del raxis.units
    raxis.id = 'member'

    # Pack
    packer = Packer(data, nordim=False, logger=logger)

    # Unpacked
    unpacked = packer.unpack(packer.packed_data)

    # Repack
    repacked = packer.repack(data)

    # Checks
    svalid = ~data[0].mask
    cdata = data.asma().reshape(data.shape[0], -1).compress(svalid.ravel(), axis=1)
    cnorm = cdata.std()
    cdata -= cdata.mean(axis=0)
    cdata /= cnorm
    assert_allclose(packer.good, svalid)
    assert_allclose(packer.sshape, data.shape[1:])
    assert_allclose(packer.norm, data.asma().std())
    assert_allclose(packer.mean, data.mean(axis=0))
    assert_allclose(packer.packed_data.shape, (svalid.sum(), data.shape[0]))
    assert_allclose(packer.packed_data, cdata.T)
    assert_allclose(unpacked, data)
    assert_allclose(repacked, packer.packed_data)

    return packer

def test_pack_mv2_scattered_without_record_fixed_norm():

    # Fake data
    lons, lats, data = create_mv2_scattered_xyzt(nt=0)
    data[:, 3:5] = MV2.masked
    data[:] *= data
    norm = data.std()*2.

    # Pack
    packer = Packer(data, nordim=True, logger=logger, mean=False, norm=norm)

    # Unpacked
    unpacked = packer.unpack(packer.packed_data)

    # Repack
    repacked = packer.repack(data)

    # Checks
    svalid = ~data.mask
    cdata = data.compressed()
    cnorm = cdata.std()*2
    cdata /= cnorm
    assert_allclose(packer.good, svalid)
    assert_allclose(packer.sshape, data.shape)
    assert_allclose(packer.norm, data.asma().std()*2)
    assert_allclose(packer.mean, 0)
    assert_allclose(packer.packed_data.shape, svalid.sum())
    assert_allclose(packer.packed_data, cdata)
    assert_allclose(unpacked, data)
    assert_allclose(repacked, packer.packed_data)

    return packer

if __name__=='__main__':
    packer = test_pack_mv2_curved_with_record()
    packer = test_pack_mv2_scattered_without_record_fixed_norm()
    print 'Done'
