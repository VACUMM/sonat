"""Commandline user interface module"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .__init__ import sonat_help, get_logger
from .config import (parse_args_cfg, get_cfg_xminmax, get_cfg_yminmax,
    get_cfg_tminmax, get_cfg_path)
from .ens import generate_pseudo_ensemble


def main():

    # Generate parser
    parser = ArgumentParser('SONAT command line interface')

    # Subparsers
    subparsers = parser.add_subparsers(title='subcommands',
        description='use "<subcommand> --help" to have more help')

    # - open help
    hparser = subparsers.add_parser('open_help', help='open the sonat help url')
    hparser.add_argument('text', help='text to search for', nargs='?')
    hparser.set_defaults(func=open_help)

    # - ensemble
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
    epparser.add_argument('ncfile', nargs='?',
        help='ensemble netcdf file')
    pgparser.set_defaults(func=ens_plot_diags)

    # Read/check config and parse commandline arguments
    args, cfg = parse_args_cfg(parser)


## HELP

def open_help(args):
    """open_help subcommand"""
    sonat_help(args.text)


## ENSEMBLE

ERRMSG_NOMODFILES = ('No model file specified. Please specify it as arguments '
    'to the command line, or in the configuration file')

def ens_gen_pseudo(parser, args):
    """ens gen_pseudo subcommand"""
    # Get the config
    cfg = args.vacumm_cfg

    # List of model files from args and config
    ncmodfiles = (args.ncmodfile if args.ncmoddile else
        get_cfg_path(cfg, 'ens', 'ncmodfiles'))
    cfg['ens']['ncmodfiles'] = ncmodfiles
    if not ncmodfiles:
        parser.error(ERRMSG_NOMODFILES)

    # Execute using config only
    ens_gen_pseudo_from_cfg(cfg)

def ens_gen_pseudo_from_cfg(cfg):
    """Take model output netcdf files and create a ensemble netcdf file"""
    # Config
    cfgd = cfg['domain']
    cfge = cfg['ens']

    # Logger
    logger = get_logger(cfg=cfg)

    # Arguments from options
    ncfile = get_cfg_path(cfg, 'ens', 'ncfile')
    ncmodfiles = get_cfg_path(cfg, 'ens', 'ncmodfiles')
    if not ncmodfiles:
        logger.error(ERRMSG_NOMODFILES)
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    time = get_cfg_tminmax(cfg, bounds=False)
    nens = cfge['nens']
    enrich = cfge['enrich']
    norms = cfge['norms']
    depths = cfge['depths']
    varnames = cfge['varnames'] or None
    getmodes = enrich > 1

    # Run and save
    ensvars = generate_pseudo_ensemble(ncmodfiles, nrens=nens, enrich=enrich,
        norms=None, lon=lon, lat=lat, time=time, depths=depths, varnames=varnames,
        getmodes=getmodes, logger=logger, asdicts=False, anomaly=True,
        ncensfile=ncfile)


def ens_plot_diags(args):
    """ens plot_diags subcommand"""
    ens_plot_diags_from_cfg(args.vacumm_cfg)

def ens_plot_diags_from_cfg(cfg):
    # Config
    cfgd = cfg['domain']
    cfge = cfg['ens']

    # Logger
    logger = get_logger('ENSPLOTDIAGS', cfg=cfg)

    # Arguments from options
    ncfile = (args.ncfile if args.ncfile else
        get_cfg_path(cfg, 'ens', 'ncfile'))
    if not ncfile:
        msg = ('No ensemble file specified. Please specify it as the first argument '
            'to the command line, or in the configuration file')
        logger.error(msg)
        parser.error(msg)
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    time = get_cfg_tminmax(cfg)
    nens = cfge['nens']
    enrich = cfge['enrich']
    norms = cfge['norms']
    depths = cfge['depths']
    varnames = cfge['varnames'] or None
    getmodes = enrich > 1
    figpatslice = get_cfg_path(cfg, 'ens', 'figpatslice')
    figpatgeneric = get_cfg_path(cfg, 'ens', 'figpatgeneric')
    htmlfile = get_cfg_path(cfg, 'ens', 'htmlfile')

    # Setup ensemble from file
    ens = Ensemble.from_file(ncfile, varnames=varnames, logger=logger)

    # Plot diags
    ens.plot_diags(htmlfile=htmlfile, figpatslice=figpatslice,
        figpatgeneric=figpatgeneric)



if __name__=='__main__':
    main()
