"""Commandline user interface module"""
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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pylab import register_cmap, get_cmap
from vcmq import dict_merge, itv_intersect

from .__init__ import sonat_help, get_logger, SONATError
from .config import (parse_args_cfg, get_cfg_xminmax, get_cfg_yminmax,
    get_cfg_tminmax, get_cfg_path)
from .misc import interpret_level
from .obs import load_obs_platform, ObsManager
from .ens import generate_pseudo_ensemble, Ensemble
from .my import load_user_code_file, SONAT_USER_CODE_FILE


def main():

    # Generate parser
    parser = ArgumentParser('SONAT command line interface')

    # Subparsers
    subparsers = parser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # Help
    hparser = subparsers.add_parser('open_help', help='open the sonat help url')
    hparser.add_argument('text', help='text to search for', nargs='?')
    hparser.set_defaults(func=open_help)


    # Ensemble
    eparser = subparsers.add_parser('ens', help='ensemble tools')
    esubparsers = eparser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # - ensemble gen_pseudo
    egparser = esubparsers.add_parser('gen_pseudo',
        help='generate a pseudo-ensemble from model outputs')
    egparser.add_argument('ncmodfile', nargs='*',
        help='model netcdf file path or pattern')
    egparser.set_defaults(func=ens_gen_pseudo)

    # - ensemble plot_diags
    epparser = esubparsers.add_parser('plot_diags',
        help='make and plot ensemble diagnostics')
    epparser.add_argument('ncensfile', nargs='?',
        help='ensemble netcdf file')
    egparser.set_defaults(func=ens_plot_diags)


    # Obs
    oparser = subparsers.add_parser('obs', help='observations tools')
    osubparsers = oparser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # - ensemble plot_diags
    opparser = osubparsers.add_parser('plot',
        help='plot observations locations or errors')
    opparser.add_argument('ncensfile', nargs='?',
        help='ensemble netcdf file')
    ogparser.set_defaults(func=obs_plot)

    # Read/check config and parse commandline arguments
    args, cfg = parse_args_cfg(parser)


## HELP

def open_help(args):
    """open_help subcommand"""
    sonat_help(args.text)


## ENSEMBLE

def ens_gen_pseudo_from_args(parser, args):
    """ens gen_pseudo subcommand"""
    # Get the config
    cfg = args.vacumm_cfg

    # List of model files from args and config
    ncmodfiles = (args.ncmodfile if args.ncmodfile else
        get_cfg_path(cfg, 'ens', 'ncmodfiles'))
    cfg['ens']['gen']['ncmodfiles'] = ncmodfiles
    if not ncmodfiles:
        parser.error('No model file specified. Please specify it as arguments '
    'to the command line, or in the configuration file')

    # Execute using config only
    return ens_gen_pseudo_from_cfg(cfg)

def ens_gen_pseudo_from_cfg(cfg):
    """Take model output netcdf files and create an ensemble netcdf file"""
    # Config
    cfgd = cfg['domain']
    cfge = cfg['ens']
    cfgeg = cfge['gen']
    cfgegf = cfgeg['fromobs']

    # Init
    logger = init_from_cfg(cfg)

    # Options from config
    ncensfile = get_cfg_path(cfg, 'ens', 'ncensfile')
    ncmodfiles = get_cfg_path(cfg, ['ens', 'gen'], 'ncmodfiles')
    if not ncmodfiles:
        raise SONATError('No model file specified')
    if not ncensfile:
        raise SONATError('No ensemble file specified')
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    time = get_cfg_tminmax(cfg, bounds=False)
    nens = cfgeg['nens']
    enrich = cfgeg['enrich']
    norms = cfg['norms']
    level = interpret_level(cfgeg['level'])
    depths = cfgeg['depths']
    varnames = cfgeg['varnames'] or None
    getmodes = enrich > 1

    # Options from obs
    if cfgef['activate']:

        # Load obs manager
        obsmanager = load_obs_from_cfg(cfg)

        # Get specs
        specs = obsmanager.get_model_specs()

        # Intervals
        margin = cfgef['margin']
        if cfgef['lon'] and specs['lon']:
            olon = specs['lon']
            if margin:
                olon = (olon[0]-margin, olon[1]+margin, olon[2])
            if lon is not None and cfgef['lon']==2:
                lon = itv_intersect(olon, lon)
            else:
                lon = olon
        if cfgef['lat'] and specs['lat']:
            olat = specs['lat']
            if margin:
                olat = (olat[0]-margin, olat[1]+margin, olat[2])
            if lat is not None and cfgef['lat']==2:
                lat = itv_intersect(olat, lat)
            else:
                lat = olat

        # Varnames
        if cfgef['varnames'] and specs['varnames']:
            if varnames is None or cfgef['varnames']==1:
                varnames = cfgef['varnames']
            else:
                varnames = list(set(varnames + specs['varnames']))

        # Depths
        olevel = interpret_level(specs['depths'])
        if cfgef['level']==1: # from obs only

            level = olevel

        elif cfgef['level']==2: # merge

            level = dict_merge(olevel, level, mergetuples=True,
                unique=True, cls=dict)
            level.setdefault('__default__', "3d") # default defaults to 3d

    # Run and save
    generate_pseudo_ensemble(ncmodfiles, nrens=nens, enrich=enrich,
        norms=None, lon=lon, lat=lat, time=time, level=level, depths=depths,
        varnames=varnames,
        getmodes=getmodes, logger=logger, asdicts=False, anomaly=True,
        ncensfile=ncensfile)

    return ncensfile


def ens_plot_diags_from_args(args):
    """ens plot_diags subcommand"""
    # Get the config
    cfg = args.vacumm_cfg

    # List of model files from args and config
    ncensfile = (args.ncensfile if args.ncensfile else
        get_cfg_path(cfg, 'ens', 'ncensfile'))
    cfg['ens']['ncensfile'] = ncensfile
    if not ncensfile:
        parser.error('No ensemble file specified. Please specify it as an argument '
    'to the command line, or in the configuration file')

    # Execute using config only
    return ens_plot_diags_from_cfg(CFG)


def ens_plot_diags_from_cfg(cfg):

    # Config
    cfgd = cfg['domain']
    cfge = cfg['ens']
    cfged = cfge['diags']
    cfgc = cfg['cmaps']
    cfgef = cfge['fromobs']
    cfgedp = cfged['plots']
#    cfgp = cfg['plots']
    cfgps = cfgp['sections']

    # Init
    logger = init_from_cfg(cfg)

    # Options from config
    ncensfile = get_cfg_path(cfg, 'ens', 'ncensfile')
    if not ncensfile:
        raise SONATError('No ensemble file specified')
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    varnames = cfge['varnames'] or None
    figpatslice = get_cfg_path(cfg, 'ens', 'figpatslice')
    figpatgeneric = get_cfg_path(cfg, 'ens', 'figpatgeneric')
    htmlfile = get_cfg_path(cfg, 'ens', 'htmlfile')
    depths = interpret_level(cfgeds('depths').dict())
    zonal_sections =  cfgps('zonal')
    merid_sections =  cfgps('merid')
    kwargs = cfged.dict().copy()
    del kwargs['plots']
    props = cfgedp.dict()
    for vname, pspecs in cfgedp.keys(): # default colormaps
        pspecs.setdefault('cmap', cfgc[vname])

    # Setup ensemble from file
    ens = Ensemble.from_file(ncensfile, varnames=varnames, logger=logger,
        lon=lon, lat=lat)

    # Plot diags
    return ens.export_diags(htmlfile, figpat_slice=figpatslice,
        figpat_generic=figpatgeneric, depths=depths, props=props,
        zonal_sections=zonal_sections, merid_sections=merid_sections,
        **kwargs)


## OBS

def load_obs_from_cfg(cfg):
    """Setup an :class:`~sonat.obs.ObsManager` using the configuration"""
    # Logger
    logger = get_logger(cfg=cfg)
    logger.verbose('Loading observations')

    # Loop on platform types
    obsplats = []
    for platform_name, platform_section in cfg['obs']['platforms'].items():

        logger.debug('Loading platform named: ' + platform_name)
        logger.debug(' Type: '+platform_section['type'])
        pfile = platform_section['file']
        logger.debug(' File: '+pfile)
        logger.debug(' Variable names: {}'.format(platform_section['varnames']))

        # Check file
        if not pfile or not os.path.exists(pfile):
            raise SONATError('Observation platform file not found: ' +
                pfile)

        # Arguments
        kwargs = platform_section.copy()
        kwargs['name'] = platform_name

        # Load
        obs = load_obs_platform(platform_section['type'], pfile, **kwargs)
        obsplats.append(obs)
    if not obsplats:
        raise SONATError('No observation platform to load were found!')

    # Init manager
    manager = ObsManager(obsplats)

    return manager


def obs_plot_from_cfg(cfg):

    # Config
    cfgo = cfg['obs']
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    cfgos = cfgo['platforms']
    cfgop = cfgo['plots']
    cfgp = cfg['plots']
    cfps = cfgp['sections']

    # Init
    logger = init_from_cfg(cfg)

    # Load obs manager
    obsmanager = load_obs_from_cfg(cfg)

    # Var names
    varnames = []
    if cfgop['locations']:
        varnames.append('locations')
    varnames.extend(cfgop['varnames'])

    # Plot
    obsmanager.plot(varnames=varnames, figpat=cfgo['figpat'],
                    lon=lon, lat=lat,
                    full3d=cfgp['full3d'], full2d=cfgp['full2d'],
                    surf=cfgp['surf'], bottom=cfgp['bottom'],
                    zonal_sections=cfgps['zonal'],
                    merid_sections=cfgps['merid'],
                    horiz_sections=cfgps['horiz'],
                    )


## MISC

def load_my_sonat_from_cfg(cfg, logger):

    # My file
    myfile = cfg['session']['usercodefile']
    myfile = myfile or SONAT_USER_CODE_FILE

    # Load it
    logger.debug('Load user code file: '+myfile)
    load_user_code_file(myfile)
    logger.verbose('Loaded user code file: '+myfile)


def register_cmaps_from_cfg(cfg, logger):
    """Register cmap aliases into matplotlib"""
    logger.debug('Registering colormaps')
    for name, cmap in cfg['cmaps']:
        register_cmap(name, cmap)


def init_from_cfg(cfg, logger):
    """Init stuff that is always performed


    Return
    ------
    logger
    """
    # Logger
    logger = get_logger(cfg=cfg)

    # Colormaps
    register_cmaps_from_cfg(cfg, logger)

    # User stuff
    load_my_sonat_from_cfg(cfg, logger)

    return logger

if __name__=='__main__':
    main()
