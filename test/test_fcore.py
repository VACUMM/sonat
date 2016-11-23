"""Test script for module :mod:`pyarm._fcore`

These tests are taken from the examples of the SANGOMA library.
"""
import os
import sys
import numpy as N

from util import THISDIR, assert_allclose

from pyarm._fcore import f_arm, f_eofcovar, f_sampleens, f_computeensstats

def test_fcore_arm():
    """Test the :func:`pyarm._fcore.f_arm` function"""

    # Sample states
    ssamples = N.array([[0.24, -0.38, 0.14], [0.51, -0.75, 0.24]], order='F')

    # Observation operator
    H = N.array([[0.0, 1.0], [1.0, 0.0]], order='F')

    # Generate data-space samples
    dsamples = N.dot(H, ssamples)

    # Observation error
    R = N.asfortranarray(N.diag([0.1, 0.1]))

    # Call ARM
    ndof = min(dsamples.shape)
    spect, U, rho_mu, status = f_arm(ndof, ssamples, dsamples, R)

    # True values
    spect_true = N.loadtxt(os.path.join(thisdir, 'sangoma_example_arm_output--RM_spectrum'),
        usecols=[1])
    U_true =  N.array([N.loadtxt(os.path.join(thisdir,
            'sangoma_example_arm_output--mode{:04d}'.format(i+1)), usecols=[1])
        for i in xrange(dsamples.shape[0])]).T
    rho_mu_true = N.array([N.loadtxt(os.path.join(thisdir,
            'sangoma_example_arm_output--modrep{:04d}'.format(i+1)), usecols=[1])
        for i in xrange(ssamples.shape[0])]).T

    # Checks
    assert_allclose(status, 0)
    assert_allclose(spect, spect_true)
    assert_allclose(U, U_true)
    assert_allclose(rho_mu, rho_mu_true)

def test_fcore_eofcovar():
    """Test the :func:`pyarm._fcore.f_eofcovar` function"""

    # Inits
    nfiles = 5
    nstate = 4
    inpath = 'inputs'
    infile = 'fieldA_{it}.txt'
    outfile_eof = 'sangoma_eofcovar_uv_eof_{}.txt'      # Files holding EOFs
    outfile_svals = 'sangoma_eofcovar_uv_svals.txt'      # Files holding singular values
    outfile_mstate = 'sangoma_eofcovar_uv_meanstate.txt' # Files holding mean state
    states = N.zeros((nstate, nfiles), order='f')
    dim_fields = N.array([1])
    offsets = 1
    remove_mstate = 1
    do_mv = 0
    meanstate = N.zeros(states.shape[0], order='f')

    # Read states
    for it in range(1, nfiles+1):
        states[:, it-1] = N.loadtxt(os.path.join(thisdir, inpath, infile).format(**locals()))

    # EOF decomposition
    stddev, svals, svec, status = f_eofcovar(dim_fields, offsets, remove_mstate,
        do_mv, states, meanstate)

    # True values
    stddev_true = 1.
    svals_true = N.loadtxt(os.path.join(thisdir, outfile_svals))
    svec_true = N.array([N.loadtxt(os.path.join(thisdir, outfile_eof.format(i+1)))
        for i in xrange(nfiles-1)]).T
    meanstate_true = N.loadtxt(os.path.join(thisdir, outfile_mstate))

    # Checks
    assert_allclose(status, 0)
    assert_allclose(stddev, stddev_true)
    assert_allclose(svals[:nfiles-1], svals_true)
    assert_allclose(svec[:, :nfiles-1], svec_true)
    assert_allclose(meanstate, meanstate_true)

def test_fcore_sampleens():
    """Test the :func:`pyarm._fcore.f_sampleens` function"""

    # Inits
    neofs = 4
    dim_state = 8
    infile_eof = 'sangoma_eofcovar_mv_eof_{}.txt'             # Files holding EOFs
    infile_svals = 'sangoma_eofcovar_mv_svals.txt'      # Files holding singular values
    infile_mstate = 'sangoma_eofcovar_mv_meanstate.txt' # Files holding mean state
    outfile_ens = 'sangoma_sampleens_ensB_{}.txt'           # Files holding ensemble states
    flag = N.array(0)

    # Read eofcovar_mv results
    svals = N.loadtxt(os.path.join(thisdir, infile_svals))
    svecs = N.array([N.loadtxt(os.path.join(thisdir, infile_eof.format(i+1)))
        for i in xrange(neofs)]).T
    meanstate = N.loadtxt(os.path.join(thisdir, infile_mstate))

    # Generate ensemble
    ens = f_sampleens(svecs, svals, meanstate, flag)

    # True values
    ens_true = N.array([N.loadtxt(os.path.join(thisdir, outfile_ens.format(i+1)))
        for i in xrange(ens.shape[1])]).T

    # Check
    assert_allclose(ens, ens_true)

def test_fcore_computeensstats():
    """Test the :func:`pyarm._fcore.f_computeens` function"""

    # Inits
    nfiles = 5
    nstate = 4
    inpath = 'inputs'
    infile = 'fieldA_{it}.txt'
    states = N.zeros((nstate, nfiles), order='f')
    element = 0

    # Read states
    for it in range(1, nfiles+1):
        states[:, it-1] = N.loadtxt(os.path.join(thisdir, inpath, infile).format(**locals()))

    # Anomaly
    meanstate = states.mean(axis=1)
    states -= meanstate.reshape(1, nstate).T

    # Compute stats
    skewness, kurtosis, status = f_computeensstats(element, meanstate, states)

    # Checks
    assert_allclose(status, 0)
    assert_allclose([skewness, kurtosis], [-1.2529E+00, -1.3728E+00], rtol=1e-4)

if __name__=='__main__':
    test_fcore_arm()
    test_fcore_eofcovar()
    test_fcore_sampleens()
    test_fcore_computeensstats()
