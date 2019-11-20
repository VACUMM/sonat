.. _scripting:

Scripting
#########

This section provides a limited view of how the
script with :mod:`sonat`.
The :ref:`unit tests <lib.test>` section also provides some examples.
For more details, have a look at the full API.

Have a look at the observation specifications
=============================================

Read the sub-sampled bathymetry for future uses

>>> import cdms2
>>> f = cdms2.open('bathy.brittany.nc')
>>> bathy = f(''elevation', lon=slice(0, 3), lat=slice(0, 3))
>>> f.close()

Load profiles from a netcdf file that contains error variables.

>>> from sonat.obs import NcObsPlatform
>>> obs = NcObsPlatform('obs.profiles.nc', name='profiles')


Plot the locations and temperature errors in 3D, at the surface
and in a zonal sections, with a marker size of 40.

>>> obs.plot(varnames=['locations', 'temp'], full3d=True, full2d=False,
... surf=True, map_margin=0.05, zonal_sections=[47.4], size=40, bathy=bathy)

.. figure:: ../../test/sonat.obs.generic_profiles_locations_map_surf.png
    :align: center

    Surface view of locations only

.. figure:: ../../test/sonat.obs.generic_profiles_temp_map_3d.png
    :align: center

    Temperature errors in 3d view

.. figure:: ../../test/sonat.obs.generic_profiles_temp_zonal_47.40.png
    :align: center

    Temperature errors in a zonal view at 47.4°N.

Print some info that are useful for model read,
like the spatial extensions and variables names.

>>> print obs.get_model_specs()
{'lat': (47.299999999999997, 48.100000000000009, 'cce'), 'depths':    id: depth
   Designated a level axis.
   units:
   Length: 21
   First:  -100.0
   Last:   0.0
   Other axis attributes:
      long_name: Levels
      axis: Z
      realtopology: linear
   Python id:  0x7f7a62792750
, 'lon': (-5.7999999999999998, -2.7999999999999998, 'cce'), 'varnames': ['temp', 'sal']}

These info can be retreive directly with some methods or
parameters like :meth:`~sonat.obs.NcObsPlatform.get_lon`,
:attr:`~sonat.obs.NcObsPlatform.depths` or
:attr:`~sonat.obs.NcObsPlatform.varnames`.

See :class:`sonat.obs.NcObsPlatform` for more
capabilities.


Generate a pseudo-ensemble from model outputs
=============================================

Generate and save the ensemble from two successive files
:file:`manga-2014-01-01.nc` and :file:`manga-2014-01-01.nc`

>>> from sonat.ens import generate_pseudo_ensemble
>>> generate_pseudo_ensemble('manga-{date:%Y-%m-%d}.nc', nrens=14, enrich=2,
...        varnames=['temp', 'sal'], time=('2014-01-01 13', '2014-01-25 12'),
...        dtfile=(15, 'day'), getmodes=True,
...        level={'temp': ('3d', 'surf')}, depths=[-40., -30, -20, -10, 0.],
...        ncensfile='ens.nc')

This function also accepts the ``lon=`` and ``lat=`` keywords to restrict
the area.
These intervals and other specifications may be retreived from observations as
explained in the previous section.

See :func:`~sonat.ens.generate_pseudo_ensemble`.


Load the ensemble
=================

Load it from a file, selecting some variable names
without their ``"_error"`` suffix.

>>> from sonat.ens import Ensemble
>>> ens = Ensemble.from_file('ens.nc', varnames=['temp', 'temp_surf', 'sal'])

or directly from variables:

>>> ens = Ensemble([temp_error, temp_surf_error, sal_error])
>>> print ens.varnames
['temp', 'temp_surf', 'sal']


Get ensemble diagnostics as a :class:`dict`

>>> diags = ens.get_diags()
>>> print diags.keys()
['mean', 'variance', 'spectrum', 'explained_variance', 'skew',
    'kurtosis', 'skewtest', 'kurtosistest', 'normaltest']
>>> diags['variance'][0].info()
*** Description of Slab temp_variance ***
id: temp_variance
shape: (4, 16, 34)
...

Save variance and spectrum to netcdf

>>> f = cdms2.open('ens.diags.nc' ,'w')
>>> for var in [diags['spectrum']] + diags['variance']:
... f.write(var)
>>> f.close()

Make some plots and get figure file names as a :class:`dict` tree

>>> figs = ens.plot_diags(variance=True, mean=False,
...        zonal_sections=[47.5], merid_sections=[-4],
...        kurtosis=False, normaltest=False, skewtest=False,
...        kurtosistest=False, skew=True)
>>> print figs['Skew']['Temp']['Map']['Surf']
sonat.ens.skew_temp_surf_map_surf.png

The tree branch types are the following:

- Diagnostic
- Variable
- Type of slice
- Location of slice

.. figure:: ../../test/sonat.ens.skew_temp_surf_map_surf.png
    :align: center

    Ensemble temperature skewness

If you pass the ``obs=obs`` argument, observation locations
will be added to the plot.

Export as html

>>> ens.export_html_diags('ens.html')


Setup an observation manager
============================

An :class:`sonat.obs.ObsManager` manages several :class:`sonat.obs.NcObsPlatform`
at the same time.

Load it from :class:`~sonat.obs.NcObsPlatform` instances

>>> obsp = NcObsPlatform('obs.profiles.nc')
>>> obsh = NcObsPlatform('obs.hfradars.nc')
>>> obss = NcObsPlatform('obs.satsst.nc')
>>> obsmanager = ObsManager([obsp, obsh, obss])

Get info

>>> print obsmanager.varnames
['v', 'u', 'temp', 'sal']
>>> print obsmanager.get_lon()
(-7.0, -2.5)
>>> obsmanager[0].name
obs.profiles

Make plots like with a single :class:`~sonat.obs.NcObsPlatform`

>>> obsmanager.plot(full2d=False, full3d=True, surf=True, bathy=bathy)

.. figure:: ../../test/sonat.obs.locations_map_surf.png
    :align: center

    Locations of observation platforms near the surface


Perform an ARM analysis
=======================

Setup the :class:`sonat.arm.ARM` instance with an
:class:`~sonat.obs.ObsManager` and an :class:`sonat.ens.Ensemble`

>>> from sonat.arm import ARM
>>> arm = ARM(ens, obsmanager, syncnorms=True, bathy=bathy)

Check input matrice to ARM analysis

>>> print arm.Af.shape # (ens_space, ens_modes)
(3902, 14)
>>> print arm.Yf.shape # (obs_space, ens_modes)
(206, 14)
>>> print arm.R.shape # (obs_space, obs_space) diagonal
(206, 206)

Analyse

>>> arm.analysis()

.. note:: The analysis is made automatically when trying to access results.

Get results as attributes :attr:`~sonat.arm.ARM.spect`, :attr:`~sonat.arm.ARM.arm`
and :attr:`~sonat.arm.ARM.rep`.

The are has many array mode variables as observation variables,
and has many modal representer variables as ensemble variables.

>>> print arm.spect.getAxisIds() # ARM spectrum
['mode']
>>> print arm.arm[0][0].getAxisIds() # Array modes dimensions of first variable of first platform
['mode', 'depth', 'station']
>>> print arm.rep[0].getAxisIds() # Modal representer of first variable
['mode', 'depth', 'lat_a', 'lon_a']

Raw matrices are also accessible with ``raw_`` attribute prefix.
These are the direct results from the :f:subr:`f_arm` subroutine,
i.e pure numeric array with no formatting and unpacking/unstacking.

>>> print arm.raw_arm.shape
(206, 14)

Plot them

>>> arm.plot_spect(shade=True, score='fnev')
>>> arm.plot_arm() # accepts same slice arguments as for obs plots
>>> arm.plot_rep(add_obs=True) # accepts same slice arguments as for ensemble plots
>>> arm.plot() # them all

.. figure:: ../../test/sonat.arm.spect.png
    :align: center

    ARM spectrum with shading below 1.

.. figure:: ../../test/sonat.arm.mode01_temp_map_surf.png
    :align: center

    First array mode for temperature at surface

.. figure:: ../../test/sonat.arm.rep.mode01_temp_surf_map_surf.png
    :align: center

    First modal representer for temperature at surface

Get a scores:

>>> arm.get_score('fnev')
4.7854420232607815
>>> arm.get_scores()
{'nev': 4, 'fnev': 4.7854420232607815, 'relvar': 67.340041766875657}

Get it from function and raw results

>>> from sonat.arm import arm_score_fnev
>>> print arm_score_fnev(arm.raw_spect, arm.raw_arm,arm.raw_rep)
4.7854420232607815


Perform a sensitivity analysis to observation locations
=======================================================

To perform such sensitivity analysis, use a :class:`sonat.arm.XYLocARMSA` instance.
This works only if a platform is mobile.

Setup the sensitivity analyser

>>> armsa = XYLocARMSA(arm)

Analyse with direct and indirect estimate

>>> resd = armsa.analyse(direct=True)
>>> resi = armsa.analyse(direct=False)

Plot results with two sort of scores

>>> armsa.plot(score_type='fnev')
>>> armsa.plot(score_type='relvar')

.. figure:: ../../test/sonat.armsa.xyloc.fnev.direct.png
    :align: center

    Sensitivity analysis with a ``fnev`` score type.

.. figure:: ../../test/sonat.armsa.xyloc.relvar.direct.png
    :align: center

    Sensitivity analysis with a ``relvar`` score type.

