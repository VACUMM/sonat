.. _ui.config.system:

The configuration system
========================

The configuration system is based on file of configuration specifications :file:`sonat.ini`.
See `configobj <http://www.voidspace.org.uk/python/configobj.html>`_ /
`validate <http://www.voidspace.org.uk/python/validate.html>`_
and :class:`vacumm.misc.config.ConfigManager`.

A :ref:`specification file <appendix.cfgspecs>`
define available options, with default values and type checks.
The user can provide a configuration file like
:ref:`this one <appendix.samples.config>`
to overwrite the default values provided by the specifications.
If a given option is not set by the user,
the default value is used, and when the option is set,
it checks against its type, and possibly other criteria.
The user file is by default :file:`sonat.cfg`,
and can be specified from commandline with
option :option:`--cfgfile <sonat --cfgfile>`.

In addition, commandline options are generated, based on
the specifications.
For instance, option :confopt:`[plots][sections]lonintervalwidth` generates
the commandline option :option:`--plots-sections-lonintervalwidth <sonat --plots-sections-lonintervalwidth>`.
When passed, the value of this option overwrites the one
from the user configuration file, which overwrites the
default configuration value.

