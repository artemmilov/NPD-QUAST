import argparse
import configparser
import os
import sys

from rdkit import RDLogger

import decoys
from nerpa_logger import NerpaLogger


# class _ToolListAction(argparse.Action):
#     def __init__(
#             self,
#             option_strings,
#             dest,
#             default=False,
#             required=False,
#             help=None,
#     ):
#         super(_ToolListAction, self).__init__(
#             option_strings=option_strings,
#             dest=dest,
#             nargs=0,
#             const=True,
#             default=default,
#             required=required,
#             help=help,
#         )
#
#     def __call__(self, parser, namespace, values, option_string=None):
#         config = configparser.ConfigParser()
#         config.read('npd_quast.ini')
#         print('\n'.join(npd_quast.SUPPORTED_TOOLS.keys()))
#         parser.exit()


def parse_args(args=None):
    """Parse arguments."""
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description="""This tool is purposed for performance 
measuring of small molecule identifiers and comparing them 
to each other. """
    )
    # parser.add_argument(
    #     '-ML',
    #     '--methods_list',
    #     action=_ToolListAction,
    #     help='show all supported tools and exit'
    # )
    parser.add_argument(
        '-M',
        '--method',
        type=str,
        help='show all supported tools and exit'
    )
    parser.add_argument(
        '-F',
        '--folder',
        type=str,
        help='working folder',
    )
    parser.add_argument(
        '-P',
        '--percent',
        type=float,
        help='percent of real spectra = count of decoys',
    )
    return parser.parse_args(args)


def handle_args(options, logger):
    mass_spectra = []
    for challenge in os.listdir(os.path.join(options.folder, 'challenges')):
        i = 0
        n = len(os.listdir(os.path.join(options.folder, 'challenges', challenge, 'spectra')))
        print('-' * 15 + 'reading mass spectra' + '-' * 15)
        for specter in os.listdir(os.path.join(options.folder, 'challenges', challenge, 'spectra')):
            if int(i / n * 49) > int((i - 1) / n * 49):
                print('#' * (int(i / n * 49) - int((i - 1) / n * 49)), end='', flush=True)
            i += 1
            mass_spectra.append(decoys.MassSpecter(file=os.path.join(
                options.folder, 'challenges', challenge, 'spectra', specter)))
        print('#\n' + '-' * 23 + 'done' + '-' * 23 + '\n')
    for challenge in os.listdir(os.path.join(options.folder, 'challenges')):
        k = int((options.percent / 100) *
                len(os.listdir(os.path.join(options.folder, 'challenges', challenge, 'spectra'))))
        print('-' * 18 + 'making decoys' + '-' * 19)
        if not os.path.exists(os.path.join(options.folder, 'challenges', challenge, 'decoys')):
            os.mkdir(os.path.join(options.folder, 'challenges', challenge, 'decoys'))
        if options.method == 'naive':
            for i, decoy in enumerate(decoys.handle_naive_method(mass_spectra, decoys_count=k)):
                decoy.write_to_file(os.path.join(options.folder, 'challenges', challenge,
                                    'decoys', 'decoy_{}.mgf'.format(i)))
        elif options.method == 'naive_extra':
            for i, decoy in enumerate(decoys.handle_naive_method_extra(mass_spectra, decoys_count=k)):
                decoy.write_to_file(os.path.join(options.folder, 'challenges', challenge,
                                    'decoys', 'decoy_{}.mgf'.format(i)))
        elif options.method == 'spectrum_based':
            for i, decoy in enumerate(decoys.handle_spectrum_based_method(mass_spectra, decoys_count=k)):
                decoy.write_to_file(os.path.join(options.folder, 'challenges', challenge,
                                    'decoys', 'decoy_{}.mgf'.format(i)))
        else:
            logger.error('Incorrect method!')
            return
    print('#\n' + '-' * 23 + 'done' + '-' * 23 + '\n')
    # logger.info('Done!')


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    logger = NerpaLogger()
    options = parse_args()
    #logger.info('Start running Decoys')
    handle_args(options, logger)


if __name__ == '__main__':
    main()
