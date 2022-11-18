import os
import shutil

import npd_quast.general
import npd_quast.report
from npd_quast.tools import SUPPORTED_TOOLS

VALID_DATA_FORMATS = ['csv']
VALID_SPECTRA_FORMATS = ['mgf']


class NPDQuastFolder:
    _folder = None

    def _check_challenges(self, logger):
        for challenge in os.listdir(
                os.path.join(self._folder, 'challenges'),
        ):
            challenge_dir = os.path.join(
                self._folder,
                'challenges',
                challenge,
            )
            if len(os.listdir(challenge_dir)) != 2:
                logger.info(
                    'Listdir length in {0} is {1} but should be 2'.format(
                        challenge_dir,
                        len(os.listdir(challenge_dir)),
                    )
                )
                return False
            if 'spectra' not in os.listdir(challenge_dir):
                logger.info(
                    'There is no \'spectra\' subdir in {} folder'.format(
                        challenge_dir
                    )
                )
                return False
            for specter in os.listdir(
                    os.path.join(challenge_dir, 'spectra'),
            ):
                if '.' not in specter:
                    logger.info('Specter {} does not have any format'.format(specter))
                    return False
                if specter.split('.')[1] not in VALID_SPECTRA_FORMATS:
                    logger.info('Specter {} is in invalid format'.format(specter))
                    return False
            files_with_dot = list(
                filter(
                    lambda f: '.' in f,
                    os.listdir(challenge_dir)
                )
            )
            if len(files_with_dot) != 1:
                logger.info('There is no database in {} dir'.format(challenge_dir))
                return False
            if files_with_dot[0].split('.')[1] not in VALID_DATA_FORMATS:
                logger.info('Database {} is in invalid format'.format(files_with_dot[0]))
                return False
        return True

    def _check_reports(self, logger):
        if not os.path.exists(os.path.join(self._folder, 'reports')):
            return True
        if len(os.listdir(
                os.path.join(self._folder, 'reports'),
        )) != 0:
            total_pages = 0
            about_pages = 0
            top_plots = 0
            quantiles_plots = 0
            down_pngs = 0
            to_right_pngs = 0
            npd_quast_pngs = 0
            for report in os.listdir(
                    os.path.join(self._folder, 'reports'),
            ):
                if os.path.isfile(
                        os.path.join(self._folder, 'reports', report),
                ):
                    if report == 'total_page.html':
                        total_pages += 1
                        continue
                    elif report == 'about_metrics_page.html':
                        about_pages += 1
                        continue
                    elif report == 'top_plot.png':
                        top_plots += 1
                        continue
                    elif report == 'quantiles_plot.png':
                        quantiles_plots += 1
                        continue
                    elif report == 'down.png':
                        down_pngs += 1
                        continue
                    elif report == 'to_right.png':
                        to_right_pngs += 1
                        continue
                    elif report == 'NPD-Quast.png':
                        npd_quast_pngs += 1
                        continue
                    else:
                        logger.info(report)
                        return False
                report_folder = os.path.join(
                    self._folder,
                    'reports',
                    report,
                )
                if set(os.listdir(report_folder)) != \
                        npd_quast.report.EMPTY_REPORT and \
                        set(os.listdir(report_folder)) != \
                        npd_quast.report.RAW_REPORT and \
                        set(os.listdir(report_folder)) != \
                        npd_quast.report.FULL_REPORT:
                    logger.info('Incorrect report consistence')
                    return False
            if {
                total_pages,
                about_pages,
                top_plots,
                quantiles_plots,
                down_pngs,
                to_right_pngs,
                npd_quast_pngs
            } not in [{0}, {1}]:
                logger.info('Reports dir is invalid')
                return False
        return True

    def _check_correctness(self, logger):
        if self._folder is None:
            logger.info('Folder is None')
            return False
        if not os.path.exists(self._folder):
            logger.info('{} does not exists'.format(self._folder))
            return False
        if (set(os.listdir(self._folder)) not in [
            {'challenges', 'reports', 'true_answers.txt'},
            {'challenges', 'reports', 'temp', 'true_answers.txt'},
            {'challenges', 'reports', 'debug', 'true_answers.txt'},
            {'challenges', 'reports', 'temp', 'debug', 'true_answers.txt'},
            {'challenges', 'true_answers.txt'},
            {'challenges', 'temp', 'true_answers.txt'},
            {'challenges', 'debug', 'true_answers.txt'},
            {'challenges', 'temp', 'debug', 'true_answers.txt'},
            {'challenges', 'reports', 'true_answers.txt', '.snakemake'},
            {'challenges', 'reports', 'temp', 'true_answers.txt', '.snakemake'},
            {'challenges', 'reports', 'debug', 'true_answers.txt', '.snakemake'},
            {'challenges', 'reports', 'temp', 'debug', 'true_answers.txt', '.snakemake'},
            {'challenges', 'true_answers.txt', '.snakemake'},
            {'challenges', 'temp', 'true_answers.txt', '.snakemake'},
            {'challenges', 'debug', 'true_answers.txt', '.snakemake'},
            {'challenges', 'temp', 'debug', 'true_answers.txt', '.snakemake'},
        ]):
            logger.info('Invalid folder format, cur format:\n{}'.format(os.listdir(self._folder)))
            return False
        return self._check_challenges(logger) and self._check_reports(logger)

    def _clean_temp(self):
        if os.path.isdir(os.path.join(self._folder, 'temp')):
            shutil.rmtree(os.path.join(self._folder, 'temp'))

    def _prepare_debug_folder(self, report, debug):
        if debug:
            if not os.path.exists(
                    os.path.join(self._folder, 'debug')
            ):
                os.mkdir(os.path.join(self._folder, 'debug'))
            if os.path.exists(
                    os.path.join(self._folder, 'debug', report)
            ):
                shutil.rmtree(os.path.join(self._folder, 'debug', report))
            os.mkdir(os.path.join(self._folder, 'debug', report))
        elif os.path.exists(
                os.path.join(self._folder, 'debug', report)
        ):
            shutil.rmtree(os.path.join(self._folder, 'debug', report))

    def __init__(self, folder, logger):
        if not os.path.exists(os.path.abspath(folder)):
            raise NotADirectoryError(
                '\"{0}\" is not a directory'.format(
                    os.path.abspath(folder),
                )
            )
        self._folder = os.path.abspath(folder)
        if (not self._check_correctness(logger)) \
                and (len(os.listdir(folder)) != 0):
            raise AttributeError(
                '\"{0}\" does not correspond to NPD-Quast format'.format(self._folder)
            )

    def make_tool_report(self, tool, report, specification, logger, debug):
        self._clean_temp()
        self._prepare_debug_folder(report, debug)
        tool.run(self._folder, report, specification, logger, debug)
        self._clean_temp()
        true_answers = npd_quast.general.parse_true_answers(
            os.path.join(
                self._folder,
                'true_answers.txt',
            ),
        )
        tool_answers_dict = {
            report: npd_quast.general.parse_tool_answers(
                os.path.join(
                    self._folder,
                    'reports',
                    report,
                    'tool_answers.txt',
                ),
            )
            for report in os.listdir(
                os.path.join(
                    self._folder,
                    'reports',
                )
            )
            if os.path.isdir(
                os.path.join(
                    self._folder,
                    'reports',
                    report,
                )
            ) and os.listdir(
                os.path.join(
                    self._folder,
                    'reports',
                    report,
                )
            ) != []
        }
        logger.info('Analyzing data...')
        npd_quast.report.write_report(
            self._folder,
            true_answers,
            tool_answers_dict,
        )
