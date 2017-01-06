#!/usr/bin/env python
# -*- coding: utf8 -*-
"""SONAT
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
__url__ = 'http://www.ifremer.fr/sonat'
__coryright__ = 'IFREMER'

import os
from warnings import warn
from vcmq import (Logger, )

import _fcore


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
    kwargs.setdefault('redirect_warnings', True)
    kwargs.setdefault('redirect_stdout', 'debug')
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
    kwargs.setdefault('level', level)
    kwargs.setdefault('logfile', file)

    if name is None:

        if sonat.LOGGER: # Use existing logger

            return sonat.LOGGER

        # Create new generic logger
        name = 'SONAT'
        sonat.LOGGER = Logger(name, **kwargs)
        return sonat.LOGGER

    elif isinstance(name, Logger):

        return name

    else:

         if os.path.exists(name):
             path = name
             name = path2logname(name)
         else:
             path = None

         if sonat.LOGGER: # Existing logger

             if sonat.LOGGER.logger.name != name: # Create a child

                 name = '{}.{}'.format(sonat.LOGGER.logger.name, name)
                 return SONATLogger(name, console=False, logfile=None)

             # Same logger to use it
             return sonat.LOGGER

         # New specific logger
         kwargs.setdefault('logfile', 'sonat.log')
         kwargs.setdefault('level', 'info')
         sonat.LOGGER = SONATLogger(name, **kwargs)
         if path is not None:
             sonat.LOGGER.debug('Running: '+path)
         return sonat.LOGGER

#: Current root :class:`~vacumm.misc.io.Logger` instance
LOGGER = None

def help(text=None, url=None):
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

#import plot # register cmaps
