import argparse
import configparser
import os
import sys

from rdkit import RDLogger

import npd_quast


class _VersionAction(argparse.Action):
    def __init__(
        self,
        option_strings,
        dest,
        default=False,
        required=False,
        help=None,
    ):
        super(_VersionAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            const=True,
            default=default,
            required=required,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        config = configparser.ConfigParser()
        config.read('npd_quast.ini')
        print('NPD-Quast {0}'.format(config['general']['version']))
        parser.exit()


class _ToolListAction(argparse.Action):
    def __init__(
        self,
        option_strings,
        dest,
        default=False,
        required=False,
        help=None,
    ):
        super(_ToolListAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            const=True,
            default=default,
            required=required,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        config = configparser.ConfigParser()
        config.read('npd_quast.ini')
        print('\n'.join(npd_quast.SUPPORTED_TOOLS.keys()))
        parser.exit()


def parse_args(args=None):
    """Parse arguments."""
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description="""This tool is purposed for performance 
measuring of small molecule identifiers and comparing them 
to each other. """
    )
    parser.add_argument(
        '-TL',
        '--tool_list',
        action=_ToolListAction,
        help='show all supported tools and exit'
    )
    parser.add_argument(
        '-V',
        '--version',
        action=_ToolListAction,
        help='show version and exit'
    )
    parser.add_argument(
        'tool',
        type=str,
        help='reporting tool name (required)',
    )
    parser.add_argument(
        'folder',
        type=str,
        help='working folder (required)',
    )
    return parser.parse_args(args)


def handle_args(options):
    config = configparser.ConfigParser()
    config.read('npd_quast.ini')
    add_tool_report(options)


def add_tool_report(options):
    if options.tool in npd_quast.tools.SUPPORTED_TOOLS.keys():
        tool = npd_quast.SUPPORTED_TOOLS[options.tool]()
        if os.path.isdir(options.folder):
            folder = npd_quast.NPDQuastFolder(options.folder)
            folder.make_tool_report(tool)
            print(
                '{0} report has been added!'.format(options.tool),
            )


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    options = parse_args()
    handle_args(options)


if __name__ == '__main__':
    main()
