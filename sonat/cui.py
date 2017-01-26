"""Commandline user interface module"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .__init__ import sonat_help, get_logger, SONATError
from .config import (parse_args_cfg, get_cfg_xminmax, get_cfg_yminmax,
    get_cfg_tminmax, get_cfg_path)
from .ens import generate_pseudo_ensemble, Ensemble


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
    epparser.add_argument('ncensfile', nargs='?',
        help='ensemble netcdf file')
    pgparser.set_defaults(func=ens_plot_diags)

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
    cfg['ens']['ncmodfiles'] = ncmodfiles
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

    # Logger
    logger = get_logger(cfg=cfg)

    # Arguments from options
    ncensfile = get_cfg_path(cfg, 'ens', 'ncensfile')
    ncmodfiles = get_cfg_path(cfg, 'ens', 'ncmodfiles')
    if not ncmodfiles:
        raise SONATError('No model file specified')
    if not ncensfile:
        raise SONATError('No ensemble file specified')
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
    generate_pseudo_ensemble(ncmodfiles, nrens=nens, enrich=enrich,
        norms=None, lon=lon, lat=lat, time=time, depths=depths, varnames=varnames,
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

    # Logger
    logger = get_logger(cfg=cfg)

    # Arguments from options
    ncensfile = get_cfg_path(cfg, 'ens', 'ncensfile')
    if not ncensfile:
        raise SONATError('No ensemble file specified')
    lon = get_cfg_xminmax(cfg)
    lat = get_cfg_yminmax(cfg)
    varnames = cfge['varnames'] or None
    figpatslice = get_cfg_path(cfg, 'ens', 'figpatslice')
    figpatgeneric = get_cfg_path(cfg, 'ens', 'figpatgeneric')
    htmlfile = get_cfg_path(cfg, 'ens', 'htmlfile')
    depths = cfg2depth(cfged.pop('depths'))
    kwargs = cfged.copy()

    # Setup ensemble from file
    ens = Ensemble.from_file(ncensfile, varnames=varnames, logger=logger,
        lon=lon, lat=lat)

    # Plot diags
    return ens.export_diags(htmlfile, figpat_slice=figpatslice,
        figpat_generic=figpatgeneric, depths=depths, **kwargs)


## MISC

def cfg2depth(depth):
    """Convert depth option to valid depth

    List are returned as lists.
    Non strings are returned whithout change.
    Strings are lower cased.
    Special "surf" and "bottom" values are returned whithout change.
    Others are converted to floats.
    """

    # From list
    if isinstance(depth, list):
        for i, d in enumerate(depth):
            depth[i] = cfg2depth(d)
        return depth

    # Scalar
    if not isinstance(depth, basestring):
        return depth
    depth = depth.lower()
    if depth not in ('bottom', 'surf'):
        depth = float(depth)
    return depth

if __name__=='__main__':
    main()
