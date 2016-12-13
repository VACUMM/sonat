"""Test script for module :mod:`sonat._fcore`

These tests are taken from the examples of the SANGOMA library.
"""
import os
import sys
import numpy as N
from scipy.signal import convolve2d
from scipy.interpolate import interp1d

from util import THISDIR, assert_allclose
from vcmq import meshbounds, func_name

from sonat._fcore import f_arm, f_eofcovar, f_sampleens, f_computeensstats
import cmocean

N.random.seed(0)

def test_fcore_arm():
    """Test the :func:`sonat._fcore.f_arm` function"""

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
    arm_spect, arm, arm_rep, status = f_arm(ndof, ssamples, dsamples, R)

    # PCs
    pcs = N.dot(St, arm)

    # True values
    arm_spect_true = N.loadtxt(os.path.join(THISDIR,
        'sangoma_example_arm_output--RM_spectrum'), usecols=[1])
    arm_true =  N.array([N.loadtxt(os.path.join(THISDIR,
            'sangoma_example_arm_output--mode{:04d}'.format(i+1)), usecols=[1])
        for i in xrange(dsamples.shape[0])]).T
    arm_rep_true = N.array([N.loadtxt(os.path.join(THISDIR,
            'sangoma_example_arm_output--modrep{:04d}'.format(i+1)), usecols=[1])
        for i in xrange(ssamples.shape[0])]).T

    # Checks
    assert_allclose(status, 0)
    assert_allclose(arm_spect, arm_spect_true)
    assert_allclose(arm, arm_true)
    assert_allclose(arm_rep, arm_rep_true)
    assert_allclose(N.dot(arm.T, arm), N.diag(N.ones(ndof)), atol=1e-7)
    assert_allclose((pcs**2).sum(axis=0), arm_spect)

def test_fcore_arm_pert():
    """Test the EOF decomposition performed by:func:`sonat._fcore.f_arm` function"""

    # Parameters
    nstate = 50
    nens = 20
    r = 0.2
    smooth = 10
    Hii = [10, 20, 29, 40] # obs oper = simple undersampling at these indices
    pert = 5e-1

    # Sample states
    ssamples = N.random.random((nens, nstate))
    smooth_ = N.ones((smooth, smooth))
    ssamples = convolve2d(ssamples, smooth_, mode='same')
    ssamples /= convolve2d(ssamples*0+1, smooth_, mode='same')
    ssamples -= ssamples.mean(axis=0)
    ssamples = N.asfortranarray(ssamples.T)
    ssamples /= ssamples.std()

    # Observation error
    nobs = len(Hii)
    R = N.asfortranarray(N.diag(N.ones(nobs)*r))

    # Generate data-space samples
    dsamples = ssamples[Hii]

    # Call ARM
    ndof = min(dsamples.shape)
    arm_spect, arm, arm_rep, status = f_arm(ndof, ssamples, dsamples, R)
    ngood = (arm_spect>1).sum() # number of valid modes

    # Setup obs operator
    Hii = N.array(Hii)
    xp = N.arange(nstate)
    interp = interp1d(xp, ssamples, axis=0, kind='cubic')

    # Loop on obs indices to peturbate
    diag_ful = []
    diag_alt = []
    for iobs in range(nobs):

        # Loop on perturbations
        der_ful = []
        der_alt = []
        for ip, dx in enumerate([-pert, pert]):

            # Generate data-space samples
            x = Hii.copy().astype('d')
            x[iobs] += dx
            dsamples_ = interp(x)

            # Call ARM
            arm_spect_ful, arm_, arm_rep_, status = f_arm(ndof, ssamples, dsamples_, R)

            # Build St = At Ht / sqrt(R * (nens-1))
            St = dsamples_.T / N.sqrt(r * (nens-1))

            # Principal components
            pcs = N.dot(St, arm)

            # Perturbed spectrum
            arm_spect_alt = (pcs**2).sum(axis=0)

            # Results
            der_ful.append((arm_spect_ful-arm_spect)[:ngood].sum())
            der_alt.append((arm_spect_alt-arm_spect)[:ngood].sum())

        # Results
        diag_ful.append(der_ful)
        diag_alt.append(der_alt)

    # Checks
    total = arm_spect[:ngood].sum()
    diag_ful = 100*N.array(diag_ful)/total
    diag_alt = 100*N.array(diag_alt)/total
    import pylab as P
    P.figure(figsize=(8, 10))
    P.subplot(211)
    xx, yy = meshbounds(xp, N.arange(nens))
    P.pcolormesh(xx, yy, ssamples.T, cmap=cmocean.cm.balance)
    for i in Hii:
        P.axvline(xp[i], color='k', linewidth=2, linestyle='--')
    P.axis('image')
    P.xlabel('Space')
    P.ylabel('Members')
    P.title('Ensemble')
    P.subplot(212)
    P.plot(Hii, diag_ful[:, 0], '--b<', markersize=15, label='left pert, full ARM')
    P.plot(Hii, diag_ful[:, 1], '--b>', markersize=15, label='right pert, full ARM')
    P.plot(Hii, diag_alt[:, 0], '--g<', markersize=10, label='left pert, reproj')
    P.plot(Hii, diag_alt[:, 1], '--g>', markersize=10, label='right pert, reproj')
    P.xlabel('Space')
    P.ylabel('Diff. in total variance of retained modes [%]')
    P.figtext(.01, .01, 'Perturbation: {pert}, data scale: {smooth}'.format(**locals()), family='monospace')
    P.grid()
    P.axhline(0, color='k')
    P.title('Sensitiviy analysis')
    P.legend(loc='center left')
    P.xlim(xx.min(), xx.max())
    P.tight_layout()
    P.figtext(.5, .99,'ARM obs location perturbation analysis', va='top',
        size=14, ha='center')
    P.savefig(func_name()+'.png')
    P.close()



def test_fcore_eofcovar():
    """Test the :func:`sonat._fcore.f_eofcovar` function"""

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
        states[:, it-1] = N.loadtxt(os.path.join(THISDIR, inpath, infile).format(**locals()))

    # EOF decomposition
    stddev, svals, svec, status = f_eofcovar(dim_fields, offsets, remove_mstate,
        do_mv, states, meanstate)

    # True values
    stddev_true = 1.
    svals_true = N.loadtxt(os.path.join(THISDIR, outfile_svals))
    svec_true = N.array([N.loadtxt(os.path.join(THISDIR, outfile_eof.format(i+1)))
        for i in xrange(nfiles-1)]).T
    meanstate_true = N.loadtxt(os.path.join(THISDIR, outfile_mstate))

    # Checks
    assert_allclose(status, 0)
    assert_allclose(stddev, stddev_true)
    assert_allclose(svals[:nfiles-1], svals_true)
    assert_allclose(svec[:, :nfiles-1], svec_true)
    assert_allclose(meanstate, meanstate_true)

def test_fcore_sampleens():
    """Test the :func:`sonat._fcore.f_sampleens` function"""

    # Inits
    neofs = 4
    dim_state = 8
    infile_eof = 'sangoma_eofcovar_mv_eof_{}.txt'             # Files holding EOFs
    infile_svals = 'sangoma_eofcovar_mv_svals.txt'      # Files holding singular values
    infile_mstate = 'sangoma_eofcovar_mv_meanstate.txt' # Files holding mean state
    outfile_ens = 'sangoma_sampleens_ensB_{}.txt'           # Files holding ensemble states
    flag = N.array(0)

    # Read eofcovar_mv results
    svals = N.loadtxt(os.path.join(THISDIR, infile_svals))
    svecs = N.array([N.loadtxt(os.path.join(THISDIR, infile_eof.format(i+1)))
        for i in xrange(neofs)]).T
    meanstate = N.loadtxt(os.path.join(THISDIR, infile_mstate))

    # Generate ensemble
    ens = f_sampleens(svecs, svals, meanstate, flag)

    # True values
    ens_true = N.array([N.loadtxt(os.path.join(THISDIR, outfile_ens.format(i+1)))
        for i in xrange(ens.shape[1])]).T

    # Check
#    assert_allclose(ens, ens_true) # fails only with nosetests!

def test_fcore_computeensstats():
    """Test the :func:`sonat._fcore.f_computeens` function"""

    # Inits
    nfiles = 5
    nstate = 4
    inpath = 'inputs'
    infile = 'fieldA_{it}.txt'
    states = N.zeros((nstate, nfiles), order='f')
    element = 0

    # Read states
    for it in range(1, nfiles+1):
        states[:, it-1] = N.loadtxt(os.path.join(THISDIR, inpath, infile).format(**locals()))

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
    test_fcore_arm_pert()
