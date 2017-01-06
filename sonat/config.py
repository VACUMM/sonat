#!/usr/bin/env python
# -*- coding: utf8 -*-
"""SONAT configuration
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

import os
from matplotlib import rcParams, rc_params_from_file
from vcmq import (ConfigManager, cfgargparse,
    adatetime, get_cmap, kwfilter)


#: Config specification file
SONAT_INIFILE = os.path.join(os.path.dirname(__file__), 'sonat.ini')

#: Default user configuration file
HYCOMVALID_DEFAULT_CFGFILE = 'sonat.cfg'

#: Matplotlib default configuration file
SONAT_DEFAULT_MATPLOTLIBRC =  os.path.join(os.path.dirname(__file__), 'matplotlibrc')

#: Matplotlib user configuration file
SONAT_USER_MATPLOTLIBRC =  'matplotlibrc'

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

def get_cfg_xminmax(cfg, bounds=None, none=True):
    """Get a ``(min,max[,bb])`` longitude interval from config"""
    return _get_domain_minmax_(cfg, 'x', -720, 720, bounds, none)

def get_cfg_yminmax(cfg, bounds=None, none=True):
    """Get a ``(min,max[,bb])`` latitude interval from config"""
    return _get_domain_minmax_(cfg, 'y', -90, 90, bounds, none)

def get_cfg_zminmax(cfg, bounds=None, none=True):
    """Get a ``(min,max[,bb])`` latitude interval from config"""
    levels = cfg['domain']['levels']
    if levels:
        zmin = min(levels)
        zmax = max(levels)
    else:
        zmin = -10000
        zmax = 0
    return _get_domain_minmax_(cfg, 'z', zmin, zmax, bounds, none)

def get_cfg_tminmax(cfg, bounds='co', none=True):
    """Get a ``(min,max[,bb])`` time interval from config"""
    return _get_domain_minmax_(cfg, 't', '1950-01-01', '2100-01-01', bounds, none)

def _get_domain_minmax_(cfg, key, defmin, defmax, bounds, none=True):
    dmin = cfg['domain'][key+'min']
    dmax = cfg['domain'][key+'max']
    if int(none)==1 and dmin is None and dmax is None:
        return
    if int(none)!=2:
        if dmin is None: dmin = defmin
        if dmax is None: dmax = defmax
    itv = (dmin, dmax)
    if bounds:
         itv += bounds,
    return itv


def load_mplrc(userfile=None):
    """Load a matplotlib or default user configuration file"""
    # Load default file first
    rcParams.update(rc_params_from_file(SONAT_DEFAULT_MATPLOTLIBRC, use_default_template=False))

    # Load user file
    userfile = str(userfile)
    if userfile=='False':
        return
    if userfile=='True':
        userfile = 'None'
    if userfile=='None':
        userfile = SONAT_USER_MATPLOTLIBRC
    if not os.path.exists(userfile):
        return
    rcParams.update(rc_params_from_file(userfile, use_default_template=False))

#load_mplrc()



