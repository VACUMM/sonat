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

import re
import os
import matplotlib
from matplotlib import rcParams, rc_params_from_file
from validate import VdtTypeError, force_list
from vcmq import (ConfigManager, cfgargparse, ArgList,
    adatetime, get_cmap, kwfilter, register_config_validator)

from .misc import interpret_level


#: Config specification file
SONAT_INIFILE = os.path.join(os.path.dirname(__file__), 'sonat.ini')

#: Default user configuration file
SONAT_DEFAULT_CFGFILE = 'sonat.cfg'

def load_cfg(cfgfile):
    """Load a configuration file"""
    return ConfigManager(SONAT_INIFILE, interpolation=False).load(cfgfile)

def parse_args_cfg(parser, cfgfilter=None):
    """Generate parse arguments,
    then return parsed arguments and configuration"""
    return cfgargparse(SONAT_INIFILE, parser, cfgfile=SONAT_DEFAULT_CFGFILE,
        interpolation=False, cfgfilter=cfgfilter)

def check_cfg_aliases(cfg, param):
    for gen_param, aliases in cfg['aliases'].items():
        for alias in aliases:
            if alias==param:
                return gen_param
    return param

def get_cfg_cmap(cfg, param, check_aliases=True):
    """Get the config colormap name for a given parameter

    Checks aliases using :func:`get_cfg_aliases`.
    """

    # Default
    default_cmap = cfg['cmaps']['default']
    if default_cmap.lower()=='none':
        default_cmap = None

    # Aliases
    if check_aliases:
        param = check_cfg_aliases(cfg, param)

    # From configuration
    cmap_name = cfg['cmaps'].get(param, default_cmap)
    if cmap_name in matplotlib.cm.cmap_d.keys():
        return cmap_name
    if cmap_name in cfg['cmaps']:
        return get_cfg_cmap(cfg, cmap_name, check_aliases=check_aliases)
    return default_cmap

def get_cfg_norms(cfg):
    """Get a dict of norms

    Checks aliases using :func:`get_cfg_aliases`.
    """
    norms = cfg['norms'].dict()
    for param in norms.keys():
        if param in cfg['aliases']:
            for alias in cfg['aliases'][param]:
                norms[alias] = norms[param]
            del norms[param]
    return norms

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

def get_cfg_plot_slice_specs(cfg, exclude=None):
    """Get kwargs that specify slices for plot functions

    Parameters
    ----------
    exclude: strings
        Exclude some slices.

    Example
    -------
    >>> print get_cfg_plot_slice_specs(cfg, exclude=["full2d", "full3d"])
    {'surf': True, 'bottom': False ....}

    Return
    ------
    dict
        With at most the following keys:
        full2d, full3d, surf, bottom,
        zonal_sections, merid_section and horiz_section.


    """
    cfgp = cfg['plots']
    kwargs = dict(full3d=cfgp['full3d'], full2d=cfgp['full2d'],
        surf=cfgp['surf'], bottom=cfgp['bottom'],
        zonal_sections=cfgp['sections']['zonal'],
        merid_sections=cfgp['sections']['merid'],
        horiz_sections=cfgp['sections']['horiz'],
    )
    if isinstance(exclude, str):
        exclude = [exclude]
    if isinstance(exclude, list):
        for exc in exclude:
            if exc in kwargs:
                del kwargs[exc]
    return kwargs

def get_cfg_path(cfg, secname, pathname, format=False, *args, **kwargs):
    """Format a relative path from the config with optional subtitutions

    This path is either absolute or relative to the "wordir" entry of the
    "session" config section.

    Examples
    --------
    >>> path = get_cfg_secpath(cfg, 'sec', 'htmlfile')
    >>> path = get_cfg_secpath(cfg, ['sec', 'subsec'], 'htmlfile')

    """
    # Get the path
    if not isinstance(secname, list):
        secname = [secname]
    sec = cfg
    for sname in secname:
        sec = sec[sname]
    paths = sec[pathname]
    if paths is None:
        return

    # Workdir for relative path
    wdir = cfg['session']['workdir']
    if not wdir:
        wdir = os.getcwd()

    # Loop on paths
    al = ArgList(paths)
    opaths = []
    for path in al.get():

        # Relative path
        if not os.path.isabs(path):
            path = os.path.abspath(os.path.join(wdir, path))

        # Format
        if format:
            path = path.format(*args, **kwargs)

        opaths.append(path)

    return al.put(opaths)


def get_cfg_obs_plot_specs(cfg, prefix=None):
    cfgo = cfg['obs']
    cfgop = cfgo['plots']
    cfgp = cfg['plots']

    kw = {'color': cfgop['colorcycle'],
          'marker': cfgop['markercycle'],
          'size': cfgop['size'],
          'legend_loc': cfgop['legendloc'],
          'map_elev': cfgp['3d']['elev'],
          'map_azim': cfgp['3d']['azim'],
          'add_bathy': cfgp['add_bathy'],
         }

    if prefix:
        for key in kw.keys():
            kw[prefix+key] = kw.pop(key)
    return kw

def is_level(value, default=None):
    """Validate a string that can be evaluated"""
    interpret_level(tuple(force_list(value, 1, 3)))
    try:
        return interpret_level(tuple(force_list(value, 1, 3)))
    except:
        raise VdtTypeError(value)

register_config_validator(level=is_level)


