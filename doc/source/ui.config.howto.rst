.. _ui.config.howto:

How to configure
================

A short example is available :ref:`here <appendix.samples.config>`.

Basic setup
-----------

You can specify a working directory other than your
current directory, with the option :confopt:`[session]workdir`.
All created files with be place by default relative to this directory.

SONAT includes a logging system.
You can set the logging level with the option :confopt:`[logger]level`,
and the logging file with :confopt:`[logger]file`.
Here are the possible levels : ``debug, verbose, info, notice, warning, error``

If you want to restrict the geographical area and time range,
set bounds in section :confsec:`[domain]`.

You can provide a user code file that will be ingested
by :func:`sonat.my.load_user_code_file`
to register user defined stuff, by setting its path
with option :confopt:`[session]usercodefile`.


Ensemble generation
-------------------

To generate a pseudo-ensemble, you must specify
the list of model files or file patterns with option :confopt:`[ens][gen]ncmodfiles`,
the requested size of your ensemble with :confopt:`[ens][gen]nens`,
and the enrichment factor with :confopt:`[ens][gen]enrich`.
The list of variable is set by the :confopt:`[ens]varnames` option,
and the output ensemble file is set by :confopt:`[ens]ncensfile` 
Variable are interpolated to fixes depths that are set by option
:confopt:`[ens][gen]depths`.
You can force the extraction of some variable at special
levels by setting options in :confopt:`[ens][gen]levels`,
like ``temp=3d, surf, bottom`` which extracts the temperature in 3D and
at the surface and at the bottom.
The normalisation coefficients may forced for each variable
in the :confsec:`[norms]` section.

It is possible to infer some options from the observation specifications
by tuning the :confsec:`[ens][gen][fromobs]` section.
For instance, the observations refer to ``"temp"`` and ``"sal"`` variable,
but you also want to include ``"u"`` and ``"v"`` in state variables.
If you set the :confopt:`[ens][gen][fromobs]varnames` to ``2``,
SONAT will merge ensemble variables and observations variables.
This can be deactivated by setting it to ``0``, and you can
force the list of state variable to those observed by settting
the option to ``1``.
The same holds for the level option :confopt:`[ens][gen][fromobs]level`.
As for the longitude and latitude bounds, the intersection are
taken when the options :confopt:`[ens][gen][fromobs]lon` or
:confopt:`[ens][gen][fromobs]lat` are set to ``2``.


Ensemble diagnostics
--------------------

Once the file :confopt:`[ens]ncensfile` is created
you can peform diagnostics that are configured in the :confsec:`[ens][diags]`
section.
The diagnostics can be switch on or off by setting them like
the skew with its option :confopt:`[ens][diags]skew`.


Observations
------------

The observation platforms are configurable as subsections
of :confsec:`[obs][platforms]`, giving their name
to the platform.
For each platform, activate it with the option
:confopt:`[obs][platforms][__many__]activate` (where ``__many__``
refer to the platform/subsection name)
and its file with :confopt:`[obs][platforms][__many__]file`.
If the platform is handled by another class than the default one
(``generic``),
set the type with :confopt:`[obs][platforms][__many__]type`.
Plots can be tuned thanks to the options of section :confsec:`[obs][plots]`.


ARM analysis
------------

The default type of score and the list of computed score types
are configurable by options :confopt:`[arm][analysis]score_type`
and :confopt:`[arm][analysis]score_types`.
The selection of modes you want to plot is set by option
:confopt:`[arm][analysis][plots]modes`, starting from ``1``,
and the list variables with :confopt:`[arm][analysis][plots]varnames`.


ARM sensitivity analyses
------------------------

Each subsection of :confsec:`[arm][sa]`
configures a sensitivity analyser, and the
only mandatory option if  :confopt:`[arm][sa][__many__]activate`,
when ``__many__`` stands for the name of the analyser.

SONAT comes with the builtin ``xyloc`` (:class:`~sonat.arm.XYLocARMSA`)
analyser configurable in section :confsec:`[arm][sa][xyloc]`.
You can set the absolute X/Y coordinates perturbation amplitude
with :confopt:`[arm][sa][xyloc]pert` in meridional degrees,
the score type with :confopt:`[arm][sa][xyloc]score_type`,
and whether the ARM analysis must be recomputed after each
perturbation or not with :confopt:`[arm][sa][xyloc]direct`.


Plots
-----

The colormaps used in plots can be defined per type of variable
in section :confsec:`[cmaps]`.

The type of plot and slice you want to perform is configurable
in section :confsec:`[plots]`: :confsec:`[plots]full3d`, :confsec:`[plots]full2d`,
:confsec:`[plots]surf` and :confsec:`[plots]bottom`.
Slices are defined in the :confsec:`[plots][sections]` section
by listing target coordinates. For instance, zonal sections
can be set by giving a list of latitudes to option
:confopt:`[plots][sections]zonal`.
Slices are computed by performing interpolations on model variable
and masking observations location that are too far from the section,
with a proximity length define by :confopt:`[plots][sections]latintervalwidth`
for zonal section, etc.



