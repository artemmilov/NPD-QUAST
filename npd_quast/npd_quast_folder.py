import os
import shutil

import npd_quast.general
import npd_quast.report
from npd_quast.tools import SUPPORTED_TOOLS

VALID_DATA_FORMATS = ['csv']
VALID_SPECTRA_FORMATS = ['mgf']


class NPDQuastFolder:
    _folder = None

    def _check_correctness(self):
        return True     # Brrrrrrrrrrrrrrrrrrrrr
        if self._folder is None:
            return False
        if not os.path.exists(self._folder):
            return False
        if (set(os.listdir(self._folder)) not in [
                {'challenges', 'reports', 'true_answers.txt'},
                {'challenges', 'reports', 'temp', 'true_answers.txt'},
        ]):
            return False
        for challenge in os.listdir(
                os.path.join(self._folder, 'challenges'),
        ):
            challenge_dir = os.path.join(
                self._folder,
                'challenges',
                challenge,
            )
            if len(os.listdir(challenge_dir)) != 2:
                return False
            if 'spectres' not in os.listdir(challenge_dir):
                return False
            for spectra in os.listdir(
                    os.path.join(challenge_dir, 'spectres'),
            ):
                if '.' not in spectra:
                    return False
                if spectra.split('.')[1] not in VALID_SPECTRA_FORMATS:
                    return False
            files_with_dot = list(
                filter(
                    lambda f: '.' in f,
                    os.listdir(challenge_dir)
                )
            )
            if len(files_with_dot) != 1:
                return False
            if files_with_dot[0].split('.')[1] not in VALID_DATA_FORMATS:
                return False
            for report in os.listdir(
                    os.path.join(self._folder, 'reports'),
            ):
                if os.path.isfile(
                    os.path.join(self._folder, 'reports', report),
                ):
                    if report != 'total_report.txt':
                        return False
                    else:
                        continue
                report_folder = os.path.join(
                    self._folder,
                    'reports',
                    report
                )
                if len(os.listdir(report_folder)) != 2:
                    return False
                if 'tool_answers.txt' not in os.listdir(report_folder):
                    return False
                if report + '.txt' not in os.listdir(report_folder):
                    return False
            return True

    def _clean_temp(self):
        if os.path.isdir(os.path.join(self._folder, 'temp')):
            shutil.rmtree(os.path.join(self._folder, 'temp'))

    def __init__(self, folder):
        if not os.path.exists(os.path.abspath(folder)):
            raise NotADirectoryError(
                '{0} is not a directory'.format(
                    os.path.abspath(folder),
                )
            )
        self._folder = os.path.abspath(folder)
        if (not self._check_correctness()) \
                and (len(os.listdir(folder)) != 0):
            raise AttributeError(
                '{0} consist something else'.format(folder)
            )

    def make_tool_report(self, tool, specification=None):
        self._clean_temp()
        tool.run(self._folder, specification)
        self._clean_temp()
        true_answers = npd_quast.general.parse_true_answers(
            os.path.join(
                self._folder,
                'true_answers.txt',
            ),
        )
        tool_answers_dict = {
            tool: npd_quast.general.parse_tool_answers(
                os.path.join(
                    self._folder,
                    'reports',
                    SUPPORTED_TOOLS[tool]().name(),
                    'tool_answers.txt',
                ),
            )
            for tool in SUPPORTED_TOOLS
            if SUPPORTED_TOOLS[tool]().name() in os.listdir(
                os.path.join(
                    self._folder,
                    'reports',
                ),
            )
        }
        npd_quast.report.write_report(
            self._folder,
            true_answers,
            tool_answers_dict,
        )
