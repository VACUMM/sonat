Core algorithms
###############

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

Such pseudo-ensemble needs to be large and made
of independent states to be representative of
all the processes at work.
In order to make the evaluation process more efficient,
the ensemble can be optimally reduced while keeping
its properties :cite:`even04` : the EOF [#eof]_ and eigenvalues
of a PCA analysis of the ensemble
are used regenerate an ensemble of independent states
and of reduced size.
We introduce an "enrichment factor" :math:`e` (:math:`> 1.0`):

.. math:: e  = \frac{n_{ens}}{n^r_{ens}}

where :math:`n^r_{ens}` is the target size of the reduced ensemble :math:`A'^r`,
and :math:`n_{ens}` is the size of the ensemble :math:`A'` which is large enough
to capture the variability.
The primes refer the anomaly.


ARM analysis
============

The Array Modes analysis :cite:`leheal09` evaluates the capacity of an observation
network potentially reduce mode errors with data assimilation,
but without acutally performing any assimilation experiment.
It takes into account model error coviariances as simlated by the
ensemble :math:`A` and observation errors :math:`R`.
Here we have removed all all :math:`r` and primes for the sake of clarity,
and we introduce
some quantities, following the formalism of :cite:`lamoal16`.
The ensemble covariance is expressed by:

.. math:: P = \frac{A A^T}{n_{ens}-1}

The scaled representer matrix is defined by:

.. math:: \chi = S S^T

where:

.. math:: S = \frac{1}{m-1}R^{-1}HA

and with :math:`H` the observation operator.
:math:`HA` is the projection of the model onto the the observations,
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


Scores
======

Scores are necessary for the quantitative evaluation of the network.
They are typically based on an analysis of the spectrum.
Here are examples of score types:

- The number of eigenvalues greater than one, which is the number of significant modes.
- The variance explained by these modes.
- The relative variance explained by these modes, which is more universal.

More score types can be easily implemented.


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

SONAT implement a few of these score types, and more can
be added by the user.

Interface to SANGOMA
====================


.. rubric:: Footnotes

.. [#enkf] Ensemble Kalman Filter
.. [#eof] Empirical Orthogonal Function
