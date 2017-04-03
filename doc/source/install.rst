Installation
############

.. highlight:: bash

With :program:`conda`
=====================

You need a version of :program:`conda`.
If not, install the Miniconda version as explained on this `page <https://conda.io/docs/install/quick.html>`_.

Then, type::

    $ conda install -c vacumm -c uvcdat -c conda-forge sonat

From sources
============

Requirements
------------

SONAT requires:

- `SANGOMA <http://www.data-assimilation.net>`_ and its dependencies for the core.
- `VACUMM <http://www.ifremer.fr/vacumm>`_ and its dependencies for the python interface.

Download
--------

Using :command:`git`::

    $ git clone https://github.com/VACUMM/sonat.git

Installation
------------

Like a standard python installation::

    $ python setup.py install [options]

   

