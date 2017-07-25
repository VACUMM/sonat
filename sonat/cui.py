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


import sys
import os
from collections import OrderedDict
from argparse import ArgumentParser

import matplotlib
from pylab import register_cmap
import cdms2
from vcmq import dict_merge, itv_intersect, checkdir

from .__init__ import sonat_help, get_logger, SONATError, get_info, SONAT_INFO
from .config import (parse_args_cfg, get_cfg_xminmax, get_cfg_yminmax,
    get_cfg_tminmax, get_cfg_path, get_cfg_plot_slice_specs,
    get_cfg_cmap, get_cfg_norms, get_cfg_obs_plot_specs, rebase_cfg_paths)
from .misc import interpret_level, dicttree_relpath
from .obs import load_obs_platform, ObsManager
from .ens import generate_pseudo_ensemble, Ensemble
from .arm import (ARM, list_arm_sensitivity_analysers,
                  get_arm_sensitivity_analyser)
from .render import render_and_export_html_template
from .my import load_user_code_file, SONAT_USER_CODE_FILE

THIS_DIR = os.path.dirname(__file__)

try:
    from sonat.test import ORDERED_MODULES as TEST_ORDERED_MODULES
    TEST_FORMAT = "sonat.test.test_{}"
except:
    test_dir = os.path.join(THIS_DIR, '..')
    sys.path.insert(0, test_dir)
    from test import ORDERED_MODULES as TEST_ORDERED_MODULES
    TEST_FORMAT = "test.test_{}"


def main(args=None):

    # Generate parser
    parser = ArgumentParser(prog='sonat',
                            description='SONAT command line interface')

    # Subparsers
    subparsers = parser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # Help
    shelp = 'open the sonat help url'
    hparser = subparsers.add_parser('help', description=shelp, help=shelp)
    hparser.add_argument('text', help='text to search for', nargs='?')
    hparser.set_defaults(func=open_help)

    # Info
    shelp = 'display info about SONAT'
    iparser = subparsers.add_parser('info', description=shelp, help=shelp)
    iparser.add_argument('key', choices=SONAT_INFO.keys(),
                         help=('a specific key info to display as one of '
                                    + ' '.join(SONAT_INFO.keys())), nargs='?')
    iparser.set_defaults(func=sonat_info)


    # Ensemble
    shelp = 'ensemble tools'
    eparser = subparsers.add_parser('ens', description=shelp, help=shelp)
    esubparsers = eparser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # - ensemble gen_pseudo
    shelp = 'generate a pseudo-ensemble from model outputs'
    egparser = esubparsers.add_parser('gen_pseudo',
                                      description=shelp, help=shelp)
    egparser.add_argument('ncmodfile', nargs='*',
        help='model netcdf file path or pattern')
    egparser.set_defaults(func=ens_gen_pseudo_from_args)

    # - ensemble plot_diags
    shelp = 'make and plot ensemble diagnostics'
    epparser = esubparsers.add_parser('plot_diags',
                                      description=shelp, help=shelp)
    epparser.add_argument('ncensfile', nargs='?',
        help='ensemble netcdf file')
    epparser.set_defaults(func=ens_plot_diags_from_args)


    # Obs
    shelp = 'observations tools'
    oparser = subparsers.add_parser('obs', description=shelp, help=shelp)
    osubparsers = oparser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # - plot
    shelp = 'plot observations locations or errors'
    opparser = osubparsers.add_parser('plot',
        description=shelp, help=shelp)
#    opparser.add_argument('platform', nargs='?',
#        help='platform name')
    opparser.set_defaults(func=obs_plot_from_args)

    # ARM
    shelp = 'ARM tools'
    aparser = subparsers.add_parser('arm', description=shelp, help=shelp)
    asubparsers = aparser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # - analysis
    shelp = 'run an ARM analysis and export results'
    aaparser = asubparsers.add_parser('analysis',
        description=shelp, help=shelp)
    aaparser.set_defaults(func=arm_analysis_from_args)

    # - sensitivity analysis
    shelp = 'run an ARM sensitivity analysis and export results'
    asparser = asubparsers.add_parser('sa',
        description=shelp, help=shelp)
    asparser.add_argument('saname', nargs='?',
        help='sensitivity analyser name, like "xyloc"')
    asparser.set_defaults(func=arm_sa_from_args)

    # GUI
    import sonat.gui
    shelp = 'Graphical User Interface'
    gparser = subparsers.add_parser('gui', description=shelp, help=shelp)
    sonat.gui.populate_argparser(gparser)
    gparser.set_defaults(func=sonat.gui.run_from_args)

    # Test
    shelp = 'launch the test suite'
    tparser = subparsers.add_parser('test', description=shelp, help=shelp)
    tparser.add_argument('name', nargs='*',
        help='name of modules to test to choose within this list: ' +
            ' '.join(TEST_ORDERED_MODULES))
    tparser.set_defaults(func=test_from_args)


    # Read/check config and parse commandline arguments
    cfg, args = parse_args_cfg(parser, args=args)
    args.func(parser, args, cfg)


## HELP

def open_help(parser, args, cfg):
    """open_help subcommand"""
    sonat_help(args.text)


## ENSEMBLE

def ens_gen_pseudo_from_args(parser, args, cfg):
    """ens gen_pseudo subcommand"""

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
    """Generate a pseudo-ensemble from model outputs

    Parameters
    ----------
    cfg: configobj/dict

    Return
    ------
    string
        Path to netcdf file
    """
    # Config
    cfge = cfg['ens']
    cfgeg = cfge['gen']
    cfgegf = cfgeg['fromobs']
    cfgegl = cfgeg['levels']

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
    level = interpret_level(cfgegl.dict())
    depths = cfgeg['depths']
    varnames = cfge['varnames'] or None
    getmodes = enrich > 1

    # Options from obs
    if cfgegf['activate']:

        # Load obs manager
        obsmanager = load_obs_from_cfg(cfg)

        # Get specs
        specs = obsmanager.get_model_specs()

        # Intervals
        margin = cfgegf['margin']
        if cfgegf['lon'] and specs['lon']:
            olon = specs['lon']
            if margin:
                olon = (olon[0]-margin, olon[1]+margin, olon[2])
            if lon is not None and cfgegf['lon']==2:
                lon = itv_intersect(olon, lon)
            else:
                lon = olon
        if cfgegf['lat'] and specs['lat']:
            olat = specs['lat']
            if margin:
                olat = (olat[0]-margin, olat[1]+margin, olat[2])
            if lat is not None and cfgegf['lat']==2:
                lat = itv_intersect(olat, lat)
            else:
                lat = olat

        # Varnames
        if cfgegf['varnames'] and specs['varnames']:
            if varnames is None or cfgegf['varnames']==1:
                varnames = cfgegf['varnames']
            else:
                varnames = list(set(varnames + specs['varnames']))

        # Depths
        olevel = interpret_level(specs['depths'])
        if cfgegf['level']==1: # from obs only

            level = olevel

        elif cfgegf['level']==2: # merge

            level = dict_merge(olevel, level, mergetuples=True,
                unique=True, cls=dict)
            level.setdefault('__default__', "3d") # default defaults to 3d

    # Run and save
    generate_pseudo_ensemble(ncmodfiles, nrens=nens, enrich=enrich,
        norms=norms, lon=lon, lat=lat, time=time, level=level, depths=depths,
        varnames=varnames,
        getmodes=getmodes, logger=logger, asdicts=False, anomaly=True,
        ncensfile=ncensfile)

    return ncensfile


def ens_plot_diags_from_args(parser, args, cfg):
    """ens plot_diags subcommand"""
    # List of model files from args and config
    ncensfile = (args.ncensfile if args.ncensfile else
        get_cfg_path(cfg, 'ens', 'ncensfile'))
    cfg['ens']['ncensfile'] = ncensfile
    if not ncensfile:
        parser.error('No ensemble file specified. Please specify it as an argument '
    'to the command line, or in the configuration file')

    # Execute using config only
    return ens_plot_diags_from_cfg(cfg)


def ens_plot_diags_from_cfg(cfg):
    """Plot ensemble diangostics

    Parameters
    ----------
    cfg: configobj/dict
    add_obs: bool
        Add observations to the plots?

    Return
    ------
    string:
        Path to HTML file
    """

    # Config
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    cfge = cfg['ens']
    cfged = cfge['diags']
    cfgedp = cfged['plots']

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
    kwargs = cfged.dict().copy()
    kwslices = get_cfg_plot_slice_specs(cfg, exclude=['full2d', 'full3d'])
    kwargs.update(kwslices)
    del kwargs['plots']
    props = cfgedp.dict()
    for param, pspecs in cfgedp.items(): # default colormaps
        if isinstance(pspecs, dict):
            pspecs.setdefault('cmap', get_cfg_cmap(cfg, param))
    norms = get_cfg_norms(cfg)
    kwobsplotspecs = get_cfg_obs_plot_specs(cfg, prefix='obs_')

    # Setup ensemble from file
    ens = Ensemble.from_file(ncensfile, varnames=varnames, logger=logger,
        lon=lon, lat=lat, norms=norms)

    # Observations
    if cfgedp['addobs']:
        kwargs.update(obs=load_obs_from_cfg(cfg), **kwobsplotspecs)



    # Plot diags
    htmlfile = ens.export_html_diags(htmlfile, figpat_slice=figpatslice,
        figpat_generic=figpatgeneric, props=props,
        lon=lon, lat=lat,
        **kwargs)

    return htmlfile


## OBS

def load_obs_from_cfg(cfg):
    """Setup an :class:`~sonat.obs.ObsManager` using the configuration"""
    # Logger
    logger = get_logger(cfg=cfg)
    logger.verbose('Loading observations')

    # Loop on platform types
    obsplats = []
    weights = []
    for platform_name, platform_section in cfg['obs']['platforms'].items():

        logger.debug('Loading platform named: ' + platform_name)
        logger.debug(' Type: '+platform_section['type'])
        pfile = platform_section['file']
        logger.debug(' File: '+pfile)
        logger.debug(' Variable names: {}'.format(platform_section['varnames']))
        logger.debug(' Weight: {}'.format(platform_section['weight']))

        # Check file
        if not pfile or not os.path.exists(pfile):
            raise SONATError('Observation platform file not found: ' +
                pfile)

        # Arguments
        kwargs = platform_section.copy()
        kwargs['name'] = platform_name
        weights.append(kwargs.pop('weight'))

        # Load
        obs = load_obs_platform(platform_section['type'], pfile, **kwargs)
        obsplats.append(obs)
    if not obsplats:
        raise SONATError('No observation platform to load were found!')

    # Norms
    norms = get_cfg_norms(cfg)

    # Init manager
    manager = ObsManager(obsplats, norms=norms, weights=weights)

    return manager

def obs_plot_from_args(parser, args, cfg):
    """obs plot subcommand"""
#    # List of model files from args and config
#    platforms = args.platform if args.platform else None

    # Execute using config only
    return obs_plot_from_cfg(cfg)#, platforms=platforms)

def obs_plot_from_cfg(cfg, platforms=None):
    """Plot observations

    Parameters
    ----------
    cfg: configobj/dict
    platforms: strings, None
        Select observation platforms by their name

    Return
    ------
    string
        Path to HTML file
    """

    # Config
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    cfgo = cfg['obs']
    cfgop = cfgo['plots']
    figpat = get_cfg_path(cfg, 'obs', 'figpat')
    kwplotspecs = get_cfg_obs_plot_specs(cfg)
    kwslices = get_cfg_plot_slice_specs(cfg)
    kwargs = {}
    kwargs.update(kwslices)
    kwargs.update(kwplotspecs)
    htmlfile = get_cfg_path(cfg, 'obs', 'htmlfile')

    # Init
    logger = init_from_cfg(cfg)

    # Load obs manager
    obsmanager = load_obs_from_cfg(cfg)

    # Bathy
    bathy = read_bathy_from_cfg(cfg, logger)
    if bathy is not None:
        obsmanager.set_bathy(bathy)

    # Var names
    varnames = []
    if cfgop['locations']:
        varnames.append('locations')
    varnames.extend(cfgop['varnames'])

    # Plot
    htmlfile = obsmanager.export_html(htmlfile,
                    varnames, figpat=figpat, lon=lon, lat=lat,
                    **kwargs)

    return htmlfile

## ARM

def arm_analysis_from_args(parser, args, cfg):
    """arm analysis subcommmand"""

    # Execute using config only
    return arm_analysis_from_cfg(cfg)

def arm_analysis_from_cfg(cfg):
    """Perform an ARM analysis

    Parameters
    ----------
    cfg: configobj/dict

    Return
    ------
    string
        Path to HTML file
    """

    # Config
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    cfga = cfg['arm']
    cfgaa = cfga['analysis']
    cfgaap = cfgaa['plots']
    cfge = cfg['ens']
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    ncensfile = get_cfg_path(cfg, 'ens', 'ncensfile')
    varnames = cfge['varnames'] or None
    norms = get_cfg_norms(cfg)
    kwargs = {}
    kwargs.update(get_cfg_plot_slice_specs(cfg, prefix='arm_'))
    kwargs.update(get_cfg_plot_slice_specs(cfg, prefix='rep_',
                                           exclude=['full2d', 'full3d']))
    kwargs.update(get_cfg_obs_plot_specs(cfg, prefix='rep_obs_'))

    # Init
    logger = init_from_cfg(cfg)

    # Load ensemble
    ens = Ensemble.from_file(ncensfile, varnames=varnames, logger=logger,
        lon=lon, lat=lat)

    # Load obs manager
    obs = load_obs_from_cfg(cfg)

    # Init ARM
    arm = ARM(ens, obs, norms=norms)

    # Bathy
    bathy = read_bathy_from_cfg(cfg, logger)
    if bathy is not None:
        arm.set_bathy(bathy)

    # Analyse and export
    htmlfile = arm.export_html(get_cfg_path(cfg, ['arm', 'analysis'], 'htmlfile'),
                    spect_figfile=get_cfg_path(cfg, ['arm', 'analysis'], 'figfile_spect'),
                    arm_figpat=get_cfg_path(cfg, ['arm', 'analysis'], 'figpat_arm'),
                    rep_figpat=get_cfg_path(cfg, ['arm', 'analysis'], 'figpat_rep'),
                    spect_score=cfgaa['score_type'],
                    score_types=cfgaa['score_types'],
                    varnames=cfgaap['varnames'],
                    modes=cfgaap['modes'],
                    lon=lon, lat=lat,
                    rep_sync_vminmax=cfgaap['rep']['syncvminmax'],
                    **kwargs
                   )

    return htmlfile

def arm_sa_from_args(parser, args, cfg):
    """arm sa subcommmand"""

    # List of sensitivity analysers
    sa_names = args.saname or None
    all_sa_names = list_arm_sensitivity_analysers()
    if sa_names is not None:
        for sa_name in sa_names:
            if sa_name not in all_sa_names:
                parser.error('Invalid sensitivity analyser name: '+ sa_name +
                             '\nPlease choose one of: '+' '.join(all_sa_names))

    # Execute using config only
    return arm_sa_from_cfg(cfg, sa_names=sa_names)

def arm_sa_from_cfg(cfg, sa_names=None):
    """Perform ARM sensitivity analyses

    Parameters
    ----------
    cfg: configobj/dict
    sa_names: strings, None
        Force the list of registered sensitivity analyserrs

    Return
    ------
    string
        Path to HTML file
    """

    # Config
    cfga = cfg['arm']
    cfgas = cfga['sa']
    cfge = cfg['ens']
    cfgo = cfg['obs']
    cfgop = cfgo['plots']
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    ncensfile = get_cfg_path(cfg, 'ens', 'ncensfile')
    varnames = cfge['varnames'] or None
    norms = get_cfg_norms(cfg)
    htmlfile = get_cfg_path(cfg, ['arm', 'sa'], 'htmlfile')

    # Init
    logger = init_from_cfg(cfg)

    # Load ensemble
    ens = Ensemble.from_file(ncensfile, varnames=varnames, logger=logger,
        lon=lon, lat=lat)

    # Load obs manager
    obs = load_obs_from_cfg(cfg)

    # Init ARM
    arm = ARM(ens, obs, norms=norms, logger=logger)

    # Bathy
    bathy = read_bathy_from_cfg(cfg, logger)
    if bathy is not None:
        arm.set_bathy(bathy)

    # SA names
    if sa_names is None:
        sa_names = list_arm_sensitivity_analysers()

    # Loop on analysers
    res = OrderedDict()
    for sa_name in sa_names:

        # Format paths
        rebase_cfg_paths(cfg, ['arm', 'sa', sa_name])

        # Config
        kwargs = cfgas[sa_name].dict()
        activate = kwargs.pop('activate')
        if not activate:
            continue

        # Setup analyser
        sa = get_arm_sensitivity_analyser(sa_name, arm, logger=logger)

        # Run and export
        res.update(sa.plot(lon=lon, lat=lat,
                                  obs_color=cfgop['colorcycle'],
                                  obs_marker=cfgop['markercycle'],
                                **kwargs))

    # Paths
    checkdir(htmlfile)
    res = dicttree_relpath(res, os.path.dirname(htmlfile))
    res = {"ARM sensitivity analyses": res}

    # Render with template
    render_and_export_html_template('dict2tree.html', htmlfile,
        title="ARM sensitivity analyses", content=res)
    arm.created(htmlfile)
    return htmlfile

### Tests

def test_from_args(parser, args, cfg):

    return test_from_cfg(cfg, args.name)


def test_from_cfg(cfg, names=None):

    # Init
    logger = init_from_cfg(cfg)

    # Get the list of valid module names
    logger.verbose('Getting the list of tests')
    if names: # check
        for name in names:
            if name not in TEST_ORDERED_MODULES:
                logger.error('Invalid test name: "'+name +
                             '". Please use one of: ' + ' '.join(TEST_ORDERED_MODULES))
            names = [name for name in TEST_ORDERED_MODULES if name in names] # reorder
            noprefix = lambda x: x if not x.startswith('test_') else x[5:]
            names = map(noprefix, names)
    else: # default
        names = TEST_ORDERED_MODULES
    logger.debug('Will test: '+' '.join(names))

    # Run nose
    logger.verbose('Performing tests')
    import nose
    successes = 0
    for i, name in enumerate(names):
        logger.debug('Testing: '+name)
        test_name = TEST_FORMAT.format(name)
        logger.debug(name)
        success = nose.run(argv=['nose', test_name])
        if success:
            logger.info('Success for '+name)
        else:
            logger.warning('Error/failure for '+name)
        successes += success
    failures = (i+1) - successes
    logger.info('Tested modules with {} success(es) and {} failure(s)'.format(
                successes, failures))


## MISC

def load_my_sonat_from_cfg(cfg, logger):

    # My file
    myfile = cfg['session']['usercodefile']

    # Load it
    if myfile is None:
        logger.debug('Will try to load default user code file')
    else:
        logger.debug('Load user code file: '+myfile)
    myfile = load_user_code_file(myfile)
    if myfile is not None:
        logger.verbose('Loaded user code file: '+myfile)


def register_cmaps_from_cfg(cfg, logger):
    """Register cmap aliases into matplotlib"""
    logger.debug('Registering colormaps')
    for name in cfg['cmaps']:
        cmap_name = get_cfg_cmap(cfg, name, check_aliases=False)
        if (name != cmap_name and
            cmap_name in matplotlib.cm.cmap_d.keys()): # add an alias
            register_cmap(name, matplotlib.cm.cmap_d[cmap_name])


def init_from_cfg(cfg):
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

def read_bathy_from_cfg(cfg, logger):
    """Return a gridded bathymetry variable that is positive on sea"""
    # Config
    cfgb = cfg['bathy']
    ncfile = get_cfg_path(cfg, 'bathy', 'ncfile')
    lon = get_cfg_xminmax(cfg, bounds='cce')
    lat = get_cfg_yminmax(cfg, bounds='cce')

    # Read
    if ncfile:

        logger.debug('Reading bathymetry: '+ncfile)
        if not os.path.exists(ncfile):
            logger.error('Bathymetry file not found: '+ncfile)
        f = cdms2.open(ncfile)
        varid = cfgb['varid']
        if not varid:
            for varid in f.listvariables():
                if len(f[varid].shape)==2:
                    break
            else:
                logger.error("Can't find a bathy variable in: "+ncfile)
        elif varid not in f.listvariables():
            logger.error('Invalid id for bathy variable: '+varid)

        kw = {}
        if lon is not None:
            kw['lon'] = lon
        if lat is not None:
            kw['lat'] = lat

        bathy = f(varid, **kw)
        f.close()

        if cfgb['samp']>1:
            bathy = bathy[::cfgb['samp'], ::cfgb['samp']]

        if cfgb['positive']:
            bathy *= -1

        return bathy

def sonat_info(parser, args, cfg):
    """sonat info subcommand"""
    print get_info(args.key)

if __name__=='__main__':
    main()
