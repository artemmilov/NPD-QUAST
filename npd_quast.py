"""This tool is purposed for performance measuring of
small molecule identifiers and comparing them to each other.
"""

import argparse
import configparser
import os
import sys

from rdkit import RDLogger

import npd_quast


def parse_args(args=None):
    """Parse arguments."""
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-T",
        "--tools",
        action='store_true',
        help="show all supported tools and exit"
    )
    parser.add_argument(
        "-V",
        "--version",
        action='store_true',
        help="show version and exit"
    )
    subparser = parser.add_subparsers(
        title='commands',
        dest="commands",
    )

    new_report = subparser.add_parser(
        name="new_report",
        help="add tool report",
    )
    new_report.add_argument(
        "tool",
        metavar="tool_name",
        type=str,
        help="reporting tool name",
    )
    new_report.add_argument(
        "folder",
        metavar="folder",
        type=str,
        help="working folder",
    )

    total = subparser.add_parser(
        name="total",
        help="add total report",
    )
    total.add_argument(
        "folder",
        metavar="folder",
        type=str,
        help="working folder",
    )

    return parser.parse_args(args)


def handle_args(options):
    config = configparser.ConfigParser()
    config.read('npd_quast.ini')
    if options.version:
        print('NPD-Quast {0}'.format(config['general']['version']))
    elif options.tools:
        print('\n'.join(npd_quast.tools.SUPPORTED_TOOLS.keys()))
    elif options.commands == 'new_report':
        add_tool_report(options)
    elif options.commands == 'total':
        add_total_report(options)


def add_tool_report(options):
    if options.tool in npd_quast.tools.SUPPORTED_TOOLS.keys():
        tool = npd_quast.tools.SUPPORTED_TOOLS[options.tool]()
        if os.path.isdir(options.folder):
            folder = npd_quast.NPDQuastFolder(options.folder)
            folder.make_tool_report(tool)
            print(
                '{0} report has been added!'.format(options.tool),
            )


def add_total_report(options):
    folder = npd_quast.NPDQuastFolder(options.folder)
    folder.make_total_report()
    print('Total report has been added!')


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    options = parse_args()
    handle_args(options)


if __name__ == '__main__':
    main()
