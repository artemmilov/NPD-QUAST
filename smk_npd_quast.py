import argparse
import configparser
import sys
import os
import json

from rdkit import RDLogger

import npd_quast
from nerpa_logger import NerpaLogger

import snakemake


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
    input_config = configparser.ConfigParser()
    input_config.read('npd_quast.ini')
    if (not os.path.isdir(input_config['dependencies']['path_to_conda'])) and \
            (not os.path.isfile(input_config['dependencies']['path_to_conda'])):
        logger.error(
            'There is no path to conda in npd_quast.ini'
        )
    if options.command == 'run_n_report':
        tool_code_name = npd_quast.tools.SUPPORTED_TOOLS[options.tool]().name()
        if (not os.path.isdir(input_config['supported tools'][tool_code_name])) and \
                (not os.path.isfile(input_config['supported tools'][tool_code_name])):
            logger.error(
                'There is no path to {} tool in npd_quast.ini'.format(options.tool)
            )
        logger.info('Running: \"run_n_report\" command...')
        with open('smk/config.yaml', 'w') as config:
            config.write("report_dir: '{}'\n".format(os.path.abspath(options.folder)))
            config.write("report_name: '{}'\n".format(options.report_name))
            config.write("run_tool: '{}'\n".format(options.tool))
            config.write("challenges: [{}]\n".format(','.join(
                os.listdir(os.path.join(options.folder, 'challenges')))))
            if options.tool == 'Dereplicator+':
                config.write("options: ['--pass-to-dereplicate \" --num_hits_to_report 10 \"', '--min-score 1']\n")
                # config.write("options: [")
                # config.write(','.join([
                #     '\'{0} {1}\''.format(
                #         str(k).replace('\'', '\\\'').replace('\"', '\\\"'),
                #         str(v).replace('\'', '\\\'').replace('\"', '\\\"')
                #     ) for k, v in json.load(open('default_configurations/{}.json'.format(
                #         tool_code_name))).items()
                #     if (v is not None) and (not isinstance(v, dict))
                # ] + [
                #     '\'{0}\''.format(str(k).replace('\'', '\\\'').replace('\"', '\\\"'))
                #     for k, v in json.load(open('default_configurations/Dereplicator_plus.json')).items()
                #     if (v is None) and (not isinstance(v, dict))
                # ]))
                # config.write("]\n")
            elif options.tool == 'Sirius':
                config.write("options: ['', '--name cur_database', '', '-c 10', '--database cur_database']\n")
            config.write('tool_dir: {}\n'.format(input_config['supported tools'][tool_code_name]))
        if options.tool == 'Dereplicator+':
            snakemake.snakemake(
                snakefile=os.path.join('smk', options.tool + '.smk'),
                workdir=options.folder,
                cores=1,
                configfiles=[os.path.join('smk', 'config.yaml')],
                # debug_dag=True,
                # debug=True,
                forceall=True,
                # unlock=True
            )
        elif options.tool == 'MAGMa+':
            snakemake.snakemake(
                snakefile=os.path.join('smk', 'MAGMa+', 'build_databases.smk'),
                workdir=options.folder,
                cores=1,
                configfiles=[os.path.join('smk', 'config.yaml')],
                # debug_dag=True,
                # debug=True,
                forceall=True,
                # unlock=True
            )
            snakemake.snakemake(
                snakefile=os.path.join('smk', 'MAGMa+', 'run.smk'),
                workdir=options.folder,
                cores=1,
                configfiles=[os.path.join('smk', 'config.yaml')],
                # debug_dag=True,
                # debug=True,
                forceall=True,
                # unlock=True,
                use_conda=True
            )
        else:
            snakemake.snakemake(
                snakefile=os.path.join('smk', 'Sirius', 'build_databases.smk'),
                workdir=options.folder,
                cores=1,
                configfiles=[os.path.join('smk', 'config.yaml')],
                # debug_dag=True,
                # debug=True,
                forceall=True,
                # unlock=True
            )
            snakemake.snakemake(
                snakefile=os.path.join('smk', 'Sirius', 'run.smk'),
                workdir=options.folder,
                cores=1,
                configfiles=[os.path.join('smk', 'config.yaml')],
                # debug_dag=True,
                # debug=True,
                forceall=True,
                # unlock=True
            )
        npd_quast.compile_reports(options, logger)
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
    options = parse_args()
    logger.info('Start running NPD-Quast')
    handle_args(options, logger)


if __name__ == '__main__':
    main()

# python smk_npd_quast.py run_n_report Dereplicator+ test
# /home/artem/Programming/bioinformatics/NPD-QUAST-test/snakemake_test/new2_unexecuted
