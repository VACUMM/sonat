"""Test script for module :mod:`sonat.pack`
"""
import os
import sys
import numpy as N
import MV2

from util import (LOGGER, assert_allclose, create_mv2_gridder_xyzt,
    create_mv2_scattered_xyzt)

from sonat.pack import Packer


def test_pack_mv2_curved_with_record():

    # Fake data
    nt = 5
    data = create_mv2_gridder_xyzt(nt=nt, rotate=30)
    data.id = 'mydata'
    data.long_name = 'My long name'
    data[:, :, 3:5, 2:4] = MV2.masked
    raxis = data.getTime()
    del raxis.axis
    del raxis.units
    raxis.id = 'member'

    # Pack
    packer = Packer(data, nordim=False, logger=LOGGER)
    packed = packer.packed_data.copy()

    # Unpacked
    unpacked = packer.unpack(packed, format=2)
    unpacked2 = packer.unpack(packed[:, :nt/2])

    # Repack
    repacked = packer.repack(data)

    # Renorm and norm back
    norm = packer.norm
    packer.set_norm(norm*2)
    renormed = packer.packed_data.copy()
    packer.set_norm(norm)
    backnormed = packer.packed_data.copy()

    # Checks
    svalid = ~data[0].mask
    cdata = data.asma().reshape(data.shape[0], -1).compress(svalid.ravel(), axis=1)
    cdata -= cdata.mean(axis=0)
    cnorm = cdata.std()
    cdata /= cnorm
    dmean = data.asma().mean(axis=0)
    assert_allclose(packer.good, svalid)
    assert_allclose(packer.sshape, data.shape[1:])
    assert_allclose(packer.mean, dmean)
    assert_allclose(packer.norm, (data.asma()-dmean).std())
    assert_allclose(packed.shape, (svalid.sum(), data.shape[0]))
    assert_allclose(packed, cdata.T)
    assert_allclose(unpacked, data)
    assert unpacked.id == 'mydata'
    assert unpacked.long_name == 'My long name'
    assert_allclose(unpacked2, unpacked[:nt/2])
    assert_allclose(repacked, packer.packed_data)
    assert_allclose(renormed, packed/2)
    assert_allclose(backnormed, packed)

    return packer

def test_pack_mv2_scattered_without_record_fixed_norm():

    # Fake data
    lons, lats, data = create_mv2_scattered_xyzt(nt=0)
    data[:, 3:5] = MV2.masked
    data[:] *= data
    norm = data.std()*2.

    # Pack
    packer = Packer(data, nordim=True, logger=LOGGER, mean=False, norm=norm)
    packed = packer.packed_data.copy()

    # Unpacked
    unpacked = packer.unpack(packer.packed_data)

    # Repack
    repacked = packer.repack(data)

    # Renorm and norm back
    norm = packer.norm
    packer.set_norm(norm*2)
    renormed = packer.packed_data.copy()
    packer.set_norm(norm)
    backnormed = packer.packed_data.copy()

    # Checks
    svalid = ~data.mask
    cdata = data.compressed()
    cnorm = cdata.std()*2
    cdata /= cnorm
    assert_allclose(packer.good, svalid)
    assert_allclose(packer.sshape, data.shape)
    assert_allclose(packer.norm, data.asma().std()*2)
    assert_allclose(packer.mean, 0)
    assert_allclose(packed.shape, svalid.sum())
    assert_allclose(packed, cdata)
    assert_allclose(unpacked, data)
    assert_allclose(repacked, packer.packed_data)
    assert_allclose(renormed, packed/2)
    assert_allclose(backnormed, packed)

    return packer

if __name__=='__main__':
    packer = test_pack_mv2_curved_with_record()
    packer = test_pack_mv2_scattered_without_record_fixed_norm()
