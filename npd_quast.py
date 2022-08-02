import argparse
import configparser
import sys
import os

from rdkit import RDLogger

import npd_quast
from nerpa_logger import NerpaLogger


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
    command = parser.add_subparsers(dest='command')
    run_n_report = command.add_parser('run_n_report')
    run_n_report.add_argument(
        'tool',
        type=str,
        help='reporting tool name',
    )
    run_n_report.add_argument(
        'report_name',
        type=str,
        help='reporting folder',
    )
    run_n_report.add_argument(
        'folder',
        type=str,
        help='working folder',
    )
    run_n_report.add_argument(
        '--config',
        type=str,
        help='Json file with configuration',
    )
    run_n_report.add_argument(
        '--debug',
        action='store_true',
        help='Run data in debug version',
    )
    compile_reports = command.add_parser('compile_reports')
    compile_reports.add_argument(
        'folder',
        type=str,
        help='working folder',
    )
    return parser.parse_args(args)


def handle_args(options, logger):
    config = configparser.ConfigParser()
    config.read('npd_quast.ini')
    if (not os.path.isdir(config['dependencies']['path_to_conda'])) and \
            (not os.path.isfile(config['dependencies']['path_to_conda'])):
        logger.error(
            'There is no path to conda in npd_quast.ini'
        )
    if options.command == 'run_n_report':
        tool_code_name = npd_quast.tools.SUPPORTED_TOOLS[options.tool]().name()
        if (not os.path.isdir(config['supported tools'][tool_code_name])) and \
                (not os.path.isfile(config['supported tools'][tool_code_name])):
            logger.error(
                'There is no path to {} tool in npd_quast.ini'.format(options.tool)
            )
        logger.info('Running: \"run_n_report\" command...')
        npd_quast.run_n_report(options, logger)
        logger.info(
            '{0} report has been reported in {1}!'.format(
                options.tool,
                options.report_name,
            )
        )
    elif options.command == 'compile_reports':
        logger.info('Running: \"compile_reports\" command...')
        npd_quast.compile_reports(options, logger)
        logger.info(
            'All reports are compiled!'
        )


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    logger = NerpaLogger()
    logger.info('Start running NPD-Quast')
    options = parse_args()
    handle_args(options, logger)


if __name__ == '__main__':
    main()
