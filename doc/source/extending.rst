.. _extending:

Extending
#########

Capabilities of SONAT can be extended in several ways, 
and this requires a bit of python scripting.

Where?
======

If you use SONAT in a scripting mode,
insert your own code
or call the :func:`sonat.my.load_user_code_file` function
before any analysis.
In the latter case, your codes are stored in one or several
user code files that can be loaded with :func:`~sonat.my.load_user_code_file`.


If you use a user interface, put all your code
and a single user code file, and define the path
in your configuration with option
:confopt:`[session]usecodefile`.

What?
=====

New generic variables
---------------------

Your netcdf variables must match
registered variables as those mentionned in :vacumm:`here <appendix.cf.html>`.
Matching is performed against the :attr:`id` (:attr:`name`) for
ensemble and observation variables, and also against the
the :attr:`standard_name` attributes  for model variables.

You can register you own variables in two ways.

- By calling :func:`vacumm.data.cf.register_cf_variable`:

    >>> register_cf_variable('ssb', standard_name='sea_surface_banana')

- By calling :func:`vacumm.data.cf.register_cf_variable_from_cfg`
  with a config file name  (:file:`cf.cfg` for instance) as argument.

  .. code-block:: ini

      [variables]

          [ssb]
              standard_name=sea_surface_banana

  Then:

      >>> register_cf_variable_from_cfg('cf.cfg')

  Alternatively, you can pass the :class:`configobj.ConfigObj`, a dictionary or
  a list of lines (like ``open('cf.cfg').readlines()``).
  For instance:

      >>> mycf = {'sbb': {'standard_name':'sea_surface_banana'} }
      >>> register_cf_variable_from_cfg(mycf)

  If you are using a user code file, your variables must be declared
  in a format compatible with :func:`~vacumm.data.cf.register_cf_variable_from_cfg`,
  with the name :attr:`vacumm_cf_variables`.


New :class:`~vacumm.data.misc.dataset.Dataset` types
----------------------------------------------------

SONAT read model outputs with function :func:`vacumm.data.DS`,
which takes a file name and a dataset type as first arguments.
The default dataset type is ``"generic"``, and other are available
like ``"mars"``.
You can register your own dataset type by deriving a
class from :class:`vacumm.data.misc.dataset.Dataset`, and registering
it with :func:`vacumm.data.register_dataset`::

    class MyModel(Dataset):
        pass

    register_dataset(MyDataset, 'mymodel')

    ds = DS('myfile.nc', 'mymodel')


A simple declaration of such class in a user code file,
and it will be registered automatically.


New observation plarform types
------------------------------

New observation platform must be derived from
:class:`sonat.obs.NcObsPlatform`.

You can overwrite the :meth:`~sonat.obs.NcObsPlatform.load`
method,
or for instance the observation operator for single variable
::meth:`~sonat.obs.NcObsPlatform.project_model` which acts
as pure interpolator::

    class MyPlatform(NcObsPlatform):

        platform_type = 'myplatform'

        def load(self, **kwargs):
            ...

        def project_model(self, var, **kwargs):
            ...

Then, register it with :func:`sonat.obs.register_obs_platform`:

    >>> register_obs_platform(MyPlatform)

If such a class is declared in the user code file,
it is automatically registered.


New ARM score functions
-----------------------

Builtin ARM score functions are listed by
:mod:`sonat.arm.ARM_SCORE_FUNCTIONS`.
They consist of short name like ``"fnev"``
and a function that starts with a **fixed prefix**
:mod:`sonat.arm.ARM_SCORE_FUNCTIONS_PREFIX`
and take the raw ARM spectrum (:attr:`~sonat.arm.ARM.raw_spect`),
array mode matrix (:attr:`~sonat.arm.ARM.raw_arm`)
and modal representer matrix (:attr:`~sonat.arm.ARM.raw_rep`) as arguments.

When scripting, you can register a new score function by calling
:func:`~sonat.arm.register_arm_score_function`.
When using a user code file, declared function whose name start
with the are :attr:`~sonat.arm.ARM_SCORE_FUNCTION_PREFIX` prefix
are automatically registered.


New ARM sensitivity analysers
-----------------------------

ARM sensitivity analysers are classes derived from
:class:`sonat.arm.ARMSA`, like :class:`sonat.arm.XYLocARMSA`.
Such class is initialised with a :class:`sonat.arm.ARM` instance.

To derive a new analyser, declare such a class and
override the :class:`sonat.arm.ARM.plot`
which must make all plots and return a :class:`dict`
of figure file names.

Then new analysers are registered with
:func:`~sonat.arm.register_arm_sensitivity_analyser`.
The are automaticall registered when declared in a user code file.


