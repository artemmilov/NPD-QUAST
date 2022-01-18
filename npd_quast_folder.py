"""
Файл, в котором и работает NPDQuast. Состоит из
папки challenges/, каждая из базы данных с молекулами
(database.smt), папки с вводными спектрами (spectres/).
Каждый файл с одним спектром. В процессе работы с
инструментом "tool_name" в папку reports/ будет
добавлен репорт с названием "report_tool_name.txt".
Также будет на время появляться файл temp/
с временными вычислениями. Он будет в процессе удалён.
"""
import os
import shutil

from report import write_report

VALID_DATA_FORMATS = ['csv', 'db']
VALID_SPECTRA_FORMATS = ['mgf', 'tree', 'txt']


class NPDQuastFolder:
    _folder = None

    def _check_correctness(self):
        if self._folder is None:
            return False
        if not os.path.exists(self._folder):
            return False
        if set(os.listdir(self._folder)) \
                != {'challenges', 'reports'}:
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
                if '.' not in report:
                    return False
                if report.split('.')[1] != 'txt':
                    return False
            return True

    def _clean_temp(self):
        shutil.rmtree(os.path.join(self._folder, 'temp'))

    def __init__(self, folder):
        if not os.path.exists(folder):
            raise NotADirectoryError(
                '{0} is not a directory'.format(folder)
            )
        self._folder = folder
        if (not self._check_correctness()) \
                and (len(os.listdir(folder)) != 0):
            raise AttributeError(
                '{0} consist something else'.format(folder)
            )

    def make_tool_report(self, tool):
        tool.run(self._folder)
        self._clean_temp()
