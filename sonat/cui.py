"""Commandline user interface"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .config import parse_args_cfg


def main():

    # Generate parser
    parser = ArgumentParser('SONAT command line interface')

    # Read/check config and parse commandline arguments
    args, cfg = parse_args_cfg(parser)


if __name__=='__main__':
    main()
