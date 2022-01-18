import glob
import os
from collections import defaultdict

from general import transform_from_to
from report import write_report


class AbstractTool:
    _spectra_format = None
    _database_format = None
    _tool_name = None

    def _init_tool(self, abs_folder):
        if not os.path.isdir(os.path.join(abs_folder, 'temp')):
            os.mkdir(os.path.join(abs_folder, 'temp'))
        if len(os.listdir(os.path.join(abs_folder, 'temp'))):
            files = glob.glob(os.path.join(abs_folder, 'temp'))
            for f in files:
                os.remove(f)
        os.mkdir(os.path.join(abs_folder, 'temp', 'spectres'))
        os.mkdir(os.path.join(abs_folder, 'temp', 'results'))
        os.mkdir(os.path.join(abs_folder, 'temp', 'tool'))

    def _write_database(self, abs_folder, database):
        if not os.path.isfile(database):
            raise ValueError('{0} is not a file'.format(database))
        database_name = os.path.split(database)[-1]
        if len(database_name.split('.')) != 2:
            raise ValueError('Indefinable format')
        cur_format = database_name.split('.')[1]
        transform_from_to(
            cur_format,
            self._database_format,
            database,
            os.path.join(abs_folder, 'temp', 'database.db'),
        )

    def _write_spectres(self, abs_folder, spectres):
        files = glob.glob(os.path.join(abs_folder, 'temp', 'spectres'))
        for f in files:
            try:
                os.remove(f)
            except IsADirectoryError:
                if os.path.split(f)[-1] != 'spectres':
                    os.rmdir(f)
                else:
                    for spectra in os.listdir(f):
                        os.remove(os.path.join(f, spectra))
        if not os.path.isdir(spectres):
            raise ValueError('{0} is not a directory'.format(spectres))
        for spectra in os.listdir(spectres):
            if len(spectra.split('.')) != 2:
                raise ValueError('Indefinable format')
            name = spectra.split('.')[0]
            cur_format = spectra.split('.')[1]
            transform_from_to(
                cur_format,
                self._spectra_format,
                os.path.join(spectres, spectra),
                os.path.join(abs_folder, 'temp', 'spectres', spectra),
            )
            res_folder = os.path.join(abs_folder, 'temp', 'results', name)
            os.mkdir(res_folder)
            with open(
                    os.path.join(res_folder, 'true_answers.txt'),
                    'w',
            ) as true_answers:
                true_answers.write(
                    'Challenge-116\t7399\n',
                )
                true_answers.write(
                    'Challenge-126\t70292\n',
                )
                true_answers.write(
                    'Challenge-140\t18829\n',
                )
                true_answers.write(
                    'Challenge-206\t82863\n',
                )

    def _run_tool(self, abs_folder, specification=None):
        pass

    def _parse_output(self, abs_folder):
        pass

    def _run_challenge(self, abs_folder, challenge, specification=None):
        database = list(
            filter(
                lambda f: len(f.split('.')) == 2,
                os.listdir(challenge),
            ),
        )[0]
        self._write_database(
            abs_folder,
            os.path.join(challenge, database),
        )
        self._write_spectres(
            abs_folder,
            os.path.join(challenge, 'spectres'),
        )
        self._run_tool(abs_folder, specification)
        self._parse_output(abs_folder)

    def _make_conclusion(self, abs_folder):
        true_answers = {}
        tool_answers = defaultdict(list)
        for result in filter(
            lambda res: os.path.isdir(
                os.path.join(abs_folder, 'temp', 'results', res),
            ),
            os.listdir(
                os.path.join(abs_folder, 'temp', 'results'),
            ),
        ):
            true_answers_data_file = os.path.join(
                abs_folder,
                'temp',
                'results',
                result,
                'true_answers.txt',
            )
            with open(true_answers_data_file) as true_answers_data:
                for true_answer in true_answers_data.read().split('\n'):
                    if true_answer != '':
                        true_answers[true_answer.split('\t')[0]] =\
                            true_answer.split('\t')[1]
            tool_answers_data_file = \
                os.path.join(abs_folder, 'temp', 'results', result, 'tool_answers.txt')
            with open(tool_answers_data_file) as tool_answers_data:
                for tool_answer in tool_answers_data.read().split('\n'):
                    if tool_answer != '':
                        tool_answers[tool_answer.split('\t')[0]].append(
                            (
                                tool_answer.split('\t')[1],
                                tool_answer.split('\t')[2],
                            ),
                        )
        write_report(
            os.path.join(abs_folder, 'reports', '{0}.txt'.format(self._tool_name)),
            true_answers,
            tool_answers,
        )

    def run(self, folder, specification=None):
        abs_folder = os.path.abspath(folder)
        self._init_tool(abs_folder)
        for challenge in os.listdir(
                os.path.join(abs_folder, 'challenges'),
        ):
            self._run_challenge(
                abs_folder,
                os.path.join(abs_folder, 'challenges', challenge),
                specification,
            )
        self._make_conclusion(abs_folder)

    def name(self):
        return self._tool_name