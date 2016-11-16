#!/usr/bin/env python
# -*- coding: utf8 -*-
"""PyARM
"""
#
# Copyright IFREMER (2016-2017)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#


__author__ = 'Stephane Raynaud, Guillaume Charria, Pierre De Mey'
__email__ = "raynaud@actimar.fr, charria@actimar.fr"
__version__ = '0.0.1'
__date__ = '2016-11-01'
__url__ = 'http://www.ifremer.fr/pyarm'
__coryright__ = 'IFREMER'

import os
from warnings import warn
from matplotlib import rcParams, rc_params_from_file
from vcmq import (ConfigManager, cfgargparse, Logger, GENERIC_VAR_NAMES,
    adatetime, get_cmap, kwfilter)

import _fcore


#: Config specification file
PYARM_INIFILE = os.path.join(os.path.dirname(__file__), 'pyarm.ini')

#: Default user configuration file
HYCOMVALID_DEFAULT_CFGFILE = 'pyarm.cfg'

#: Matplotlib default configuration file
PYARM_DEFAULT_MATPLOTLIBRC =  os.path.join(os.path.dirname(__file__), 'matplotlibrc')

#: Matplotlib user configuration file
PYARM_USER_MATPLOTLIBRC =  'matplotlibrc'


class PyARMError(Exception):
    pass

class PyARMWarning(UserWarning):
    pass

def pyarm_warn(message, stacklevel=2):
    """Issue a :class:`PyARMWarning`"""
    warn(message, PyARMWarning, stacklevel=stacklevel)


class PyARMLogger(Logger):
    def created(self, msg):
        msg = 'Created: '+msg
        self.info(msg)

def get_logger(name=None, cfg=None, **kwargs):
    """Automatically setup a logger for the current script"""
    kwargs.setdefault('redirect_warnings', True)
    kwargs.setdefault('redirect_stdout', 'debug')
    import pyarm
    if cfg is not None:
        kwargs.setdefault('level', cfg['logger']['level'])
        kwargs.setdefault('logfile', cfg['logger']['file'])
    if name is None:

        if pyarm.LOGGER: # Use existing logger

            return pyarm.LOGGER

        # Create new generic logger
        name = 'PYARM'
        pyarm.LOGGER = Logger(name, **kwargs)
        return pyarm.LOGGER

    elif isinstance(name, Logger):

        return name

    else:

         if os.path.exists(name):
             path = name
             name = path2logname(name)
         else:
             path = None

         if pyarm.LOGGER: # Existing logger

             if pyarm.LOGGER.logger.name != name: # Create a child

                 name = '{}.{}'.format(pyarm.LOGGER.logger.name, name)
                 return PyARMLogger(name, console=False, logfile=None)

             # Same logger to use it
             return pyarm.LOGGER

         # New specific logger
         kwargs.setdefault('logfile', 'pyarm.log')
         kwargs.setdefault('level', 'info')
         pyarm.LOGGER = PyARMLogger(name, **kwargs)
         if path is not None:
             pyarm.LOGGER.debug('Running: '+path)
         return pyarm.LOGGER

#: Current root :class:`~vacumm.misc.io.Logger` instance
LOGGER = None

def help(text=None, url=None):
    """Open pyarm website in a web browser and optionally search for a string

    :Params:

        - **text**, optional: String to search for.
        - **recent**, optional: Use the most recent version of the documentation.
    """
    if url is None: url = __url__
    from webbrowser import open
    if text is not None:
        if not isinstance(text, basestring):
            if hasattr(text, 'func_name'):
                text = text.func_name
            else:
                text = text.__class__.__name__
        if not text.startswith('/'):
            text = '/search.html?q=%s&check_keywords=yes&area=default'%text
        url += text
    open(url,  new=2)

def get_cfg_cmap(cfg, param):
    """Get the config colormap for a given parameter"""

    # Default
    default_cmap = cfg['cmaps']['default']
    if default_cmap.lower()=='none':
        default_cmap = None

    # From configuration
    cmap_name = cfg['cmaps'].get(param, default_cmap)
    try:
        return get_cmap(default_cmap)
    except:
        return get_cmap()


def load_mplrc(userfile=None):
    """Load a matplotlib or default user configuration file"""
    # Load default file first
    rcParams.update(rc_params_from_file(PYARM_DEFAULT_MATPLOTLIBRC, use_default_template=False))

    # Load user file
    userfile = str(userfile)
    if userfile=='False':
        return
    if userfile=='True':
        userfile = 'None'
    if userfile=='None':
        userfile = PYARM_USER_MATPLOTLIBRC
    if not os.path.exists(userfile):
        return
    rcParams.update(rc_params_from_file(userfile, use_default_template=False))

#load_mplrc()


class _Base_(object):

    def __init__(self, logger=None, **kwargs):

        self.logger = get_logger(logger, **kwfilter(kwargs, 'logger_'))

        self.debug('Instantiate '+self.__class__.__name__)

import plot # register cmaps
