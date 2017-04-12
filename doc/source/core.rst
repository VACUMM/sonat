.. _core:

The SONAT core
##############


Pseudo-ensemble generation
==========================

A classic dynamical ensemble evolves with time
and is created by performing different runs of the
same model, that differ for instance by their initial
or boundary conditions, by their forcing, or by the value
of their internal parameters.
Such ensemble simulate the model errors and is used
by stochastic assimilation schemes such as the Enkf
(Ensemble Kalman Filter) :cite:`even94`.

A pseudo-ensemble is static and generated from
different model states :cite:`even03`.
It is representative of the model variance and may used
as a proxy for model error thanks to a scaling factor
between 0 an 1.
In our case, this factor is set to 1 by default,
which means that the observation network will be evaluated
against its capacity to sample the variability, not
to reduce model errors though data assimilation.

In the current implementation, model states
is loaded at regular dates
by :func:`sonat.ens.load_model_at_regular_dates`
without any time interpolation.
This ensures that states are independent from one
one another, and makes the read faster.

Such pseudo-ensemble needs to be large and made
of independent states to be representative of
all the processes at work.
In order to make the evaluation process more efficient,
the ensemble can be **optimally reduced** while keeping
its properties :cite:`even04`: the EOF [#eof]_ and eigenvalues
of a PCA [#pca]_ analysis of the ensemble
are used regenerate an ensemble of independent states
and of reduced size.
We introduce an **enrichment factor** :math:`e` (:math:`> 1.0`):

.. math:: e  = \frac{n_{ens}}{n^r_{ens}}

where :math:`n^r_{ens}` is the target size of the reduced ensemble :math:`A'^r`,
and :math:`n_{ens}` is the size of the ensemble :math:`A'` which is large enough
to capture the variability and
where the primes refer the anomaly.

The pseudo-ensemble is generated with :func:`sonat.ens.generate_pseudo_ensemble`,
to which your typically a netcdf file pattern, the size of your requested ensemble
and an enrichment factor which defaults to 1.

.. figure:: ../../test/TEST_CUI/ENS/sonat.ens.spectrum.png

    Relative spectrum of retained mode for ensemble reduction.

.. figure:: ../../test/TEST_CUI/ENS/sonat.ens.explained_variance_temp_map_surf.png

    Explained variance of temperature after ensemble reduction.



Ensemble diagnostics
====================

A few diagnostics ar available (:meth:`sonat.ens.Ensemble.get_diags`):
the ensemble variance and mean (~zero), the :func:`~scipy.stats.skew`,
the :func:`~scipy.stats.skewtest`,
the :func:`~scipy.stats.kurtosis`,
the :func:`~scipy.stats.kurtosistest` and the :func:`~scipy.stats.normaltest`.

Some other stats are related to the generation of the pseudo-ensemble:
the relative spectrum and the local variance.

.. figure:: ../../test/TEST_CUI/ENS/sonat.ens.variance_v_map_surf.png

    Ensemble variance of meridional velocity.
    
.. figure:: ../../test/TEST_CUI/ENS/sonat.ens.skew_v_map_surf.png

    Skewness of the ensemble of meridional velocity.


ARM analysis
============

The Array Modes analysis :cite:`leheal09` evaluates the capacity of an observation
network potentially reduce mode errors with data assimilation,
but without acutally performing any assimilation experiment.
It takes into account model error coviariances as simlated by the
ensemble state anomalies :math:`A` and observation error covariance matrix :math:`R`.
Here we have removed all all :math:`r` and primes for the sake of clarity,
and we introduce
some quantities, following the formalism of :cite:`lamoal16`.
The ensemble covariance is expressed by:

.. math:: P = \frac{A A^T}{n_{ens}-1}

The scaled representer matrix is defined by:

.. math:: \chi = S S^T

where:

.. math:: S = \frac{1}{\sqrt{m-1}}R^{-1}HA

is the scaled ensemble state anomalies projected onto observations,
with :math:`H` the observation operator.
:math:`HA` is the projection of the ensemble anomalies onto the the observations,
which may be in some cases more than a simple interpolation.
The matrix :math:`\chi` can be seen as measure of the covariances
relative to the observation errors.
In the case of a pseudo-ensemble, it is a signal-to-noise matrix.

The ARM analysis is based on an EOF decomposition of :math:`\chi`,
which is actually achieved through an SVD analysis of the matrix :math:`S`.

.. math:: \chi = \mu \sigma \mu^T

The observation network is quantitatively evaluated by analysis
the **spectrum** of the decomposition, especially the
of eigen values :math:`\sigma` greatee than 1.
The spatial properties the network are given by the EOF
of the decompositions, also called the **array modes** :math:`\mu`.
And the signature of these modes in the model space
are given by the **modal reprensenters**:

.. math:: \rho =  \frac{1}{m-1} A S^T \mu

The model representer show how the observational network
impact the state variables, whether they are observed or not.

The ARM analysis is performed by the :meth:`sonat.ARM.analysis` method,
which store the spectrum (:attr:`~sonat.ARM.spect`),
the array modes (:attr:`sonat.ARM.arm`) and the
modal representers (:attr:`sonat.ARM.rep`)
as formatted arrays.
Raw results and other matrices are
also available.

.. figure:: ../../test/TEST_CUI/ARM/sonat.arm.spectrum.png

    Example of ARM spectrum.
    The shaded area marks eigen values lower than 1.


ARM scores
==========

Scores are necessary for the quantitative evaluation of the network.
They are typically based on an analysis of the spectrum.
Here are examples of score types:

- The number of eigenvalues greater than one,
  which is the number of significant modes (see :func:`~sonat.arm.arm_score_nev`).
- The fractional version of the latest score type (see :func:`~sonat.arm.arm_score_fnev`).
- The variance explained by these modes (see :func:`~sonat.arm.arm_score_relvar`).
- The relative variance explained by these modes, which is more universal.

The current list of score type is accessible with
:func:`sonat.arm.list_arm_score_types`.
More score types can be easily implemented and registered
with :func:`sonat.arm.register_arm_score_function`.


Sensitivity analyses
====================

Sensitivity analyses are useful for instance to check the stability
of your score with respect to parameters, or to have a indication
of how to change your network in order to optimise it.

Indeed, two very different network can be clearly compared to assess which
one is the best.
In reality, an observational network is generally already existing,
and the goal is to setup an extension to this network with
a limited number of freedom.
You can easily test the sensitivity of the score to slight changes
in observation positions in all directions.
This may tell you whether you must change your configuration or not,
and how to do it.
The :class:`sonat.arm.XYLocARMSA` class is a sensitivity analyser
that tests the effect of infinitesimal changes in the position
of observations () that are mobile

.. figure:: ../../test/arm.sa.xyloc.fnev.png

    Sensitivity analysis to observation locations as performed
    by :class:`sonat.arm.XYLocARMSA` sensitivity analyser.

New sensitivity analysers can be implemented by inheriting
from :class:`sonat.arm.ARMSA` and registered by
:func:`sonat.arm.register_arm_sensitivity_analyser`.
The list of sensitivity analyser names is available
with :func:`sonat.arm.list_arm_sensitivity_analysers`.


Interface to SANGOMA
====================

This python function calls the :f:func:`f_arm` fortran subroutine,
which is wrapper to the SANGOMA


.. rubric:: Footnotes

.. [#enkf] Ensemble Kalman Filter
.. [#eof] Empirical Orthogonal Function
.. [#pca] Principal Component Analysis
.. [#svd] Singular Value Decomposition
