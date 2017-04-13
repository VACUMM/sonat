Conventions
###########


.. _conv.vars:

Name of variables
=================

In order to track variable between model and observation data file,
all variables must have a known generic short name, without underscores.
This helps for instance to normalise all occurences of temperature data
with the same factor, by simply refering to the "temp" variable.

The list of builtin names is available :vacumm:`here <appendix.cf.html>`.

New variables can be registered by calling
:func:`~vacumm.data.cf.register_cf_variable` or
:func:`~vacumm.data.cf.register_cf_variables_from_cfg`.


.. _conv.depths:

The special depths ``"surf"`` and ``"bottom"``
==============================================

These two special depths are used with impact on slicing and naming.
The ``"surf"`` depth refere to the top of a 3D volumic data,
and ``"bottom"`` to its bottom.
The ``"surf"`` depth may be sligthly different from an interpolation at the zero depth.
Colocating variables at the bottom requires the an estimate of the bathymetry.

Ensemble variable names may be suffixed with ``"_surf"`` or ``"_bottom"`` with
when an extraction at this depths is required.
Otherwise, variables are interpolated to a fixed number of vertical levels.
These depths are also interpreted into the observations files.
Finally, they are ones of the possibles slices used for plots.


.. _conv.mod:

Model files
===========

All files must be by default in netcdf file format and CF compliant.

.. _conv.mode.mod:

Model outputs for pseudo-ensemble generation
--------------------------------------------

Model output must be CF compliant netcdf files with known variables (see :ref:`conv.vars`).
To be recongnised, these variable must either have a known name (id) or a known standard_name attribute.
Once read by :func:`vacumm.data.DS`, it is formatted with a the matching known id, like "temp" for temperature.
Three dimensional variable are interpolated to a fixed vertical axis before other processing,
like the ensemble formation and reduction.

.. _conv.mode.ens:

Ensemble file
-------------

A netcdf ensemble file must contain all members,
and the member dimension must be the first one for all variables.
All variables must have a known name, optionally suffixed with ``"_surf"`` or ``"_bottom"``.
The vertical dimension must be in Z coordinates, positive up.

Here is an :ref:`example <appendix.samples.ens>`.

.. _conv.obs:

Observation files
=================

Observation netcdf files must contain known variables suffixed with "_error",
and all variables must have the same dimensions.

The horizontal dimensions must the last.
In case of a gridded dataset like HF radar currents or satellite data,
the last two dimensions must refer to the (Y,X) 1D axes
(:ref:`hfradars <appendix.samples.obs.hfradars>`
and :ref:`satsst <appendix.samples.obs.satsst>`).
If not gridded, the dataset contains scattered observations and
the spatial axis is interpreted as a spatial index
(:ref:`profiles <appendix.samples.obs.profiles>`).
In this case, a longitude and a latitude variables must be found.

The depths may be specified in several ways:

- A depth axis: there is vertical dimension in addition the horizontal one(s) (:ref:`profiles <appendix.samples.obs.profiles>`).
  This axis is the same to all locations, and if data are not available at given depth, they must be marked as masked.
- A depth variable: data are scattered both horizontally and vertically (like a scanfish or Seasoar).
- A "depth" file attribute: a scalar value in meters or special string value ``"surf"`` or ``"bottom"``.


A mobility specification (see :ref:`core.sa`) may be included to tell if the XY locations are
mobile for adjustment of the network or not.
It can take several forms:

- An integer variable with the same horizontal dimensions as error variables.
- An integer file attribute that applies to all locations.

An integer value of 1 makes the location mobile.
By default, a platform is not mobile.

Here are :ref:`examples of observation files <appendix.samples.obs>`.


Slicing data for plots
======================

The same formalism is used in the library or the user interfaces
to perform slices for plotting observation and ensemble data.
Model variable are sliced with :func:`sonat.misc.slice_gridded_var`
and observation variables are sliced with
:func:`sonat.misc.mask_scattered_locs`.
The following slice specifications are supported.

``horiz_section``
    It is a single or a list of negative depth floats.
    Model variables are interpolated to these depths,
    and only observation variables are within a depth interval around
    the specified depths are plotted.

``surf``
    This slice plot the top level of 3D variables, which is the case
    for instance of most model variables.
    For observation variables, if they have a ``"surf"`` depth, they are
    directly plotted, otherwise this is equivalent to ``horiz_section=0.`` .

``bottom``
    For model variables, the deepest unmasked data is kept (:func:`sonat.misc.sclice_bottom`).
    For observation variables, if they have a ``"bottom"`` depth, they are
    directly plotted, otherwise the bathymetry is required at observation locations,
    and only observations that are within a depth interval around this variable depth
    are plotted.

``zonal_sections``
    It is a single or a list of latitudes.
    Model variables are interpolated at this latitude, which
    may results in a 1D  or 2D plot.
    Observation locations that are within a laitude interval around
    the reference latitudes are plotted.

``merid_sections``
    See ``zonal_sections``, but with longitudes.
   
    
Main plots are performed by :func:`sonat.plot.plot_gridded_var` for model variable
and :func:`sonat.plot.plot_scattered_locs`.

