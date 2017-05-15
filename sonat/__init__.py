#!/usr/bin/env python
# -*- coding: utf8 -*-
"""SONAT

The following environment variables are declared

.. envvar:: SONAT_LIB_DIR

    Directory path of the library

.. envvar:: SONAT_TEST_DIR

    Directory path of the library unit tests scripts

.. envvar:: SONAT_DATA_DIR

    Directory path of the library data samples

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
__version__ = '1.0.0'
__date__ = '2017-05-15'
__url__ = 'http://relay.actimar.fr/~raynaud/sonat'
__copyright__ = 'IFREMER'

import os
from collections import OrderedDict
from warnings import warn, filterwarnings
import cdms2
from vcmq import (Logger, netcdf4)

import _fcore

cdms2.setAutoBounds('off')
netcdf4()

filterwarnings("ignore", "The get_axis_bgcolor function was deprecated")
filterwarnings("ignore", "setting an item on a masked array which has a shared mask")
filterwarnings("ignore", "axes.hold is deprecated")
filterwarnings("ignore", "The ishold function was deprecated in version 2.0")

SONAT_LIB_DIR = os.path.dirname(__file__)


#: Bottom generic variable names
BOTTOM_VARNAMES = ['turb']

class SONATError(Exception):
    pass

class SONATWarning(UserWarning):
    pass

def sonat_warn(message, stacklevel=2):
    """Issue a :class:`SONATWarning`"""
    warn(message, SONATWarning, stacklevel=stacklevel)


class SONATLogger(Logger):
    def created(self, msg):
        msg = 'Created: '+msg
        self.info(msg)

def get_logger(name=None, cfg=None, **kwargs):
    """Automatically setup a logger for the current script"""
    import sonat
    if cfg is not None:
        level = cfg['logger']['level']
        file = cfg['logger']['file']
        redirect_warnings = cfg['logger']['redirect_warnings']
        redirect_stdout = cfg['logger']['redirect_stdout']
        if str(redirect_stdout).lower()=='false':
            redirect_stdout = False
    else:
        level = 'info'
        file = None
        redirect_warnings = True
        redirect_stdout = 'debug'
    kwargs.setdefault('redirect_warnings', redirect_warnings)
    kwargs.setdefault('redirect_stdout', redirect_stdout)
    kwargs.setdefault('level', level)
    kwargs.setdefault('logfile', file)

    if name is None:

        if sonat.LOGGER: # Use existing logger

            return sonat.LOGGER

        # Create new generic logger
        name = 'SONAT'
        sonat.LOGGER = SONATLogger(name, **kwargs)
        return sonat.LOGGER

    elif isinstance(name, Logger):

        return name

    else:

        # Existing logger
        if sonat.LOGGER: 

            if sonat.LOGGER.name != name: # Create a child

                name = '{}.{}'.format(sonat.LOGGER.name, name)
                return SONATLogger(name, console=False, logfile=None)

            # Same logger to use it
            return sonat.LOGGER

        # New specific logger
        kwargs.setdefault('logfile', 'sonat.log')
        kwargs.setdefault('level', 'info')
        sonat.LOGGER = SONATLogger(name, **kwargs)
        return sonat.LOGGER

#: Current root :class:`~vacumm.misc.io.Logger` instance
LOGGER = None

def sonat_help(text=None, url=None):
    """Open sonat website in a web browser and optionally search for a string

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

def get_data_dir():
    """Get the directory that contains sample data"""
    data_dir = os.path.join(SONAT_LIB_DIR, 'data')
    if not os.path.exists(data_dir):
        data_dir = os.path.abspath(os.path.join(SONAT_LIB_DIR,  '..', 'data'))
    return data_dir

def get_test_dir():
    """Get the directory that contains test scripts"""
    test_dir = os.path.join(SONAT_LIB_DIR, 'test')
    if not os.path.exists(test_dir):
        test_dir = os.path.abspath(os.path.join(SONAT_LIB_DIR,  '..', 'test'))
    return test_dir


SONAT_INFO = OrderedDict(version=__version__, author=__author__, email=__email__,
                       date=__date__, url=__url__, copyright=__copyright__,
                       data_dir=get_data_dir(), test_dir=get_test_dir(),
                       lib_dir=SONAT_LIB_DIR)

def info(key=None):
    """Print SONAT info"""
    print get_info(key)

def get_info(key=None):
    if key:
        if key not in SONAT_INFO:
            raise SONATError('Invalid info key: {}. '.format(key) +
                             'Please choose one of: '+' '.join(SONAT_INFO.keys()))
        return SONAT_INFO[key]
    return """SONAT-{version}
  Date: {date}
  Author: {author}
  Email: {email}
  URL: {url}
  Copyright: {copyright}
  Library dir: {lib_dir}
  Data dir: {data_dir}
  Test dir: {test_dir}""".format(**SONAT_INFO)

os.environ['SONAT_LIB_DIR'] = SONAT_LIB_DIR
os.environ['SONAT_TEST_DIR'] = SONAT_INFO['test_dir']
os.environ['SONAT_DATA_DIR'] = SONAT_INFO['data_dir']
