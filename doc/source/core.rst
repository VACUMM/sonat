.. _core:

Core
####


Pseudo-ensemble generation
==========================

A classic dynamical ensemble evolves with time
and is created by performing different runs of the
same model, that differ for instance by their initial
and/or boundary conditions, by their forcing, or by the value
of their internal parameters.
Such ensemble simulate the model errors and is used
by stochastic assimilation schemes such as the Enkf [#enkf]_
:cite:`even94`.

A pseudo-ensemble is static and generated from
different model states :cite:`even03`.
It is representative of the model variance and may be used
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
This ensures that states are independent from
one another, and makes the read faster.

Such pseudo-ensemble needs to be large and made
of independent states to be representative of
all the processes at work.
In order to make the evaluation process more efficient,
the ensemble can be **optimally reduced** while keeping
its properties almost intact :cite:`even04`: the EOF [#eof]_ and eigenvalues
of a PCA [#pca]_ analysis of the ensemble
are used to regenerate an ensemble of independent states
and of reduced size.
We introduce an **enrichment factor** :math:`e` (:math:`> 1.0`):

.. math:: e  = \frac{n_{ens}}{n^r_{ens}}

where :math:`n^r_{ens}` is the target size of the reduced ensemble :math:`A'^r`,
and :math:`n_{ens}` is the size of the ensemble :math:`A'` which is large enough
to capture the variability and
where the primes refer to the anomaly.

The pseudo-ensemble is generated with :func:`sonat.ens.generate_pseudo_ensemble`,
to which you provide at least a netcdf file pattern, the size of your requested ensemble
and an enrichment factor which defaults to 1.
The core part of sample creation is performed is two main steps by fortran
subroutines which are direct wrapper to the SANGOMA library's equivalent subroutines:

- the EOF decomposition of the covariance matrix is made by :f:func:`f_eofcovar`,
  which also provides a spectrum (see :numref:`fig_core_pseudo_spect`);
- the sample generation is made by :f:func:`f_sampleens`.

Note that the generated sample has lower variance than its original field
variance since a reduced number of EOFs is kept (:numref:`fig_core_pseudo_expl`).

.. _fig_core_pseudo_spect:

.. figure:: ../../test/CUI/ENS/sonat.ens.spectrum.png
    :align: center

    Example of relative spectrum of retained mode for ensemble reduction.

.. _fig_core_pseudo_expl:

.. figure:: ../../test/CUI/ENS/sonat.ens.explained_variance_temp_map_surf.png
    :align: center

    Example of explained variance of temperature after ensemble reduction.



Ensemble diagnostics
====================

A few diagnostics ar available (:meth:`sonat.ens.Ensemble.get_diags`):
the ensemble variance and mean (~zero), the :func:`~scipy.stats.skew`,
the :func:`~scipy.stats.skewtest`,
the :func:`~scipy.stats.kurtosis`,
the :func:`~scipy.stats.kurtosistest` and the :func:`~scipy.stats.normaltest`.

Some other stats are related to the generation of the pseudo-ensemble:
the relative spectrum and the local variance.

.. figure:: ../../test/CUI/ENS/sonat.ens.variance_v_surf_map_surf.png
    :align: center

    Ensemble variance of meridional velocity.

.. figure:: ../../test/CUI/ENS/sonat.ens.variance_sal_merid_4.55W.png
    :align: center

    Ensemble variance of salinity at 4.55°W.

.. figure:: ../../test/CUI/ENS/sonat.ens.skew_u_surf_map_surf.png
    :align: center

    Skewness of the ensemble of zonal velocity.


Observation platforms
=====================

An observation platform is defined by the observations errors
of a series of variable with the same coordinates.
Locations may be randomly distributed like for profiles, ferryboxes or scanfishes,
or gridded like for satellite or HF radar data.

A plaform may also integrate an :ref:`observation operator <core.obsoper>`
that is more complex than a simple interpolation, for some variables.

.. figure:: ../../test/sonat.obs.locations_map_3d.png
    :align: center

    Example of 3D view of all platforms.

Since pseudo-ensemble are used in the ARM analysis,
a **platform-specific weight** :math:`W` may be set to take into account
the time sampling of processes.
A low weight must be defined for a platform that has
a measurement time step :math:`\tau_{obs}`
significantly greater than the processes
it is designed to observe :math:`\tau_{proc}`.
Conversely, this weight saturates (to ``1``) when the time step becomes
lower than the process time scale.

.. math:: W = max\left( \frac{\tau_{proc}}{\tau_{obs}}, 1\right)



ARM analysis
============

The Array Modes analysis :cite:`leheal09` evaluates the capacity of an observation
network potentially to reduce mode errors with data assimilation,
but without acutally performing any assimilation experiment.
It takes into account model error coviariances as simulated by the
ensemble state anomalies :math:`A` and observation error covariance matrix :math:`R`.
Here we have removed all :math:`r` and primes for the sake of clarity,
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
:math:`HA` is the projection of the ensemble anomalies onto the observations,
which may be in some cases more than a simple interpolation.
The matrix :math:`\chi` can be seen as measure of the covariances
relative to the observation errors.
In the case of a pseudo-ensemble, it is a signal-to-noise matrix.

The ARM analysis is based on an EOF decomposition of :math:`\chi`,
which is actually achieved through an SVD [#svd]_ analysis of the matrix :math:`S`.

.. math:: \chi = \mu \sigma \mu^T

The observation network is quantitatively evaluated by analysis
the **spectrum** of the decomposition (Fig. :numref:`fig_core_arm_spect`), especially the
of eigen values :math:`\sigma` greater than 1.
The spatial properties the network are given by the EOF
of the decompositions, also called the **array modes** :math:`\mu`
(Fig. :numref:`fig_core_arm_arm_temp3d` and :numref:`fig_core_arm_arm_usurf`).
And the signature of these modes in the model space
are given by the **modal representers** (Fig. :numref:`fig_core_arm_rep_temp`):

.. math:: \rho =  \frac{1}{m-1} A S^T \mu

The model representer shows how the observational network
impact the state variables, whether they are observed or not.

The ARM analysis is performed by the :meth:`sonat.ARM.analysis` method,
which store the spectrum (:attr:`~sonat.ARM.spect`),
the array modes (:attr:`sonat.ARM.arm`) and the
modal representers (:attr:`sonat.ARM.rep`)
as formatted arrays.
Raw results and other matrices are
also available.

.. _fig_core_arm_spect:

.. figure:: ../../test/CUI/ARM/sonat.arm.spect.png
    :align: center

    Example of ARM spectrum.
    The shaded area marks eigen values lower than 1.

.. _fig_core_arm_arm_temp3d:

.. figure:: ../../test/CUI/ARM/sonat.arm.arm.mode01_temp_map_3d.png
    :align: center

    Example of a 3D view of the first array mode for temperature.

.. _fig_core_arm_arm_usurf:

.. figure:: ../../test/CUI/ARM/sonat.arm.arm.mode01_u_map_surf.png

    Example of a surface view of the first array mode for zonal velocity.


.. _fig_core_arm_rep_temp:

.. figure:: ../../test/CUI/ARM/sonat.arm.rep.mode01_temp_map_surf.png
    :align: center

    Example of a surface view of the first modal representer of temperature,
    with an indication of surface observation locations.


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


.. _core.sa:

Sensitivity analyses
====================

Sensitivity analyses are useful for instance to check the stability
of your score with respect to parameters, or to have a indication
of how to change your network in order to optimise it.

Indeed, two very different networks can be clearly compared to assess which
one is the best.
In reality, an observational network is generally already existing,
and the goal is to setup an extension to this network with
a limited number of freedom.
You can easily test the sensitivity of the score to slight changes
in observation positions in all directions.
This may tell you whether you must change your configuration or not,
and how to do it.

The :class:`sonat.arm.XYLocARMSA` class is a sensitivity analyser
that tests the effect of infinitesimal changes in the
horizontal position of observations that are mobile.
A platform can be fully, partially or not mobile, and is not mobile by default.
This is specified by a ``mobility`` variable or global attribute
in an observation file (see :ref:`conv.obs`).

A ``score_type`` must be chosen, and tests may be performed
either with a ``direct`` or ``indirect`` method:
in the first case, a full ARM analysis is performed
after each perturbation, while in the latter case,
the original state variables are interpolated to
the new positions (:math:`H^*`) to give a new ensemble state
anomalies matrix :math:`H^*A`
at observation locations,
which is converted to a perturbed :math:`S^*` matrix,
which is then projected onto original array modes
(EOFs) to provide perturbed expansion coefficients :math:`a^*`:

.. math:: a^* =  S^{*T} \mu

from which a perturbed spectrum :math:`\sigma^*` is computed


.. math:: \sigma^* = \frac{1}{n_{ens}} a^{*T}a^*

This spectrum can now be used to compute scores.
The advantage of this approach is its low cost since
no SVD decomposition is needed at each perturbation.

.. _fig_core_armsa_direct:

.. figure:: ../../test/sonat.armsa.xyloc.fnev.indirect.png
    :align: center

    Sensitivity analysis to observation locations as performed
    by :class:`sonat.arm.XYLocARMSA` sensitivity analyser,
    with a ``fnev`` score type and a ``indirect`` method.

.. _fig_core_armsa_indirect:

.. figure:: ../../test/sonat.armsa.xyloc.relvar.indirect.png
    :align: center

    Same as previous figure, but with an ``relvar`` score type.

In the examples of figures :numref:`fig_core_armsa_direct`
and :numref:`fig_core_armsa_indirect`, the direct and
indirect methods show similar results.
They suggest to move the most eastern profile to the
south.

New sensitivity analysers can be implemented by inheriting
from :class:`sonat.arm.ARMSA` and registered by
:func:`sonat.arm.register_arm_sensitivity_analyser`.
The list of sensitivity analyser names is available
with :func:`sonat.arm.list_arm_sensitivity_analysers`.


Multivariate analyses
=====================

SONAT support multivariate analyses at all levels.
It is made possible thanks to normalisation factors
applied to ensemble states and observation errors.

When the ensemble states normalisation factor is not provided,
it is computed with using the standard deviation.
As for the observation errors, the normalisation
factor :math:`\sigma`
is approximated using function :func:`sonat.misc.sqrt_errors_norm`
applied to :math:`r = \sqrt{R}`.

.. math:: \sigma = \sqrt{\frac{1}{N}\sum r^{2}_{i}}

Note that, normalisation factors can provided per variable type,
like temperature or salinity.


.. _core.obsoper:

The observation operator
========================

The observation operator is in charge of the projection of
model outputs to observation locations.
It is only made of interpolation tasks in the current version of SONAT
(:meth:`sonat.obs.NcObsPlatform.project_model`).
Latter, builtin and user-made complex operators will be possible,
like for satellite SST [#sst]_, HF radar radial currents or
ocean color.


.. rubric:: Footnotes

.. [#enkf] Ensemble Kalman Filter
.. [#eof] Empirical Orthogonal Function
.. [#pca] Principal Component Analysis
.. [#sst] Sea Surface Temperature
.. [#svd] Singular Value Decomposition
