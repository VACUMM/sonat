"""Commandline user interface module"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .__init__ import sonat_help
from .config import parse_args_cfg


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

    # Read/check config and parse commandline arguments
    args, cfg = parse_args_cfg(parser)


def open_help(args):
    sonat_help(args.text)



if __name__=='__main__':
    main()
