import configparser
import os
import shutil

from npd_quast import general
from npd_quast.report import write_report


class AbstractTool:
    _spectra_format = None
    _database_format = None
    _tool_name = None
    _location = None

    def __init__(self):
        config = configparser.ConfigParser()
        config.read('npd_quast.ini')
        self._location = config['supported tools'][self._tool_name]

    def _init_tool(self, abs_folder):
        if not os.path.isdir(os.path.join(abs_folder, 'temp')):
            os.mkdir(os.path.join(abs_folder, 'temp'))
        if len(os.listdir(os.path.join(abs_folder, 'temp'))) > 0:
            for subdir in os.listdir(
                    os.path.join(abs_folder, 'temp'),
            ):
                shutil.rmtree(subdir)
        os.mkdir(os.path.join(abs_folder, 'temp', 'spectres'))
        os.mkdir(os.path.join(abs_folder, 'temp', 'tool'))
        if self._tool_name in os.listdir(
                os.path.join(abs_folder, 'reports'),
        ):
            shutil.rmtree(
                os.path.join(abs_folder, 'reports', self._tool_name),
            )
        os.mkdir(
            os.path.join(abs_folder, 'reports', self._tool_name),
        )
        os.mknod(
            os.path.join(
                abs_folder,
                'reports',
                self._tool_name,
                'tool_answers.txt',
            ),
        )

    def _convert_database(self, from_database, to_database):
        pass

    def _write_database(self, abs_folder, database):
        if not os.path.isfile(database):
            raise ValueError('{0} is not a file'.format(database))
        database_name = os.path.split(database)[-1]
        if len(database_name.split('.')) != 2:
            raise ValueError('Indefinable format')
        self._convert_database(
            database,
            os.path.join(
                abs_folder,
                'temp',
                'database.{0}'.format(
                    self._database_format,
                ),
            ),
        )

    def _convert_spectra(self, from_spectra, to_spectra):
        pass

    def _write_spectres(self, abs_folder, spectres):
        if os.path.isdir(os.path.join(abs_folder, 'temp', 'spectres')):
            shutil.rmtree(os.path.join(abs_folder, 'temp', 'spectres'))
            os.mkdir(os.path.join(abs_folder, 'temp', 'spectres'))
        if not os.path.isdir(spectres):
            raise ValueError('{0} is not a directory'.format(spectres))
        for spectra in os.listdir(spectres):
            if len(spectra.split('.')) != 2:
                raise ValueError('Indefinable format')
            self._convert_spectra(
                os.path.join(spectres, spectra),
                os.path.join(abs_folder, 'temp', 'spectres', spectra),
            )

    def _run_tool(self, abs_folder, specification=None):
        pass

    def _parse_output(self, abs_folder, challenge_name):
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
        self._parse_output(abs_folder, os.path.split(challenge)[-1])

    def _make_conclusion(self, abs_folder, default_ranks=None):
        true_answers = general.parse_true_answers(
            os.path.join(
                abs_folder,
                'true_answers.txt',
            ),
        )
        tool_answers = general.parse_tool_answers(
            os.path.join(
                abs_folder,
                'reports',
                self._tool_name,
                'tool_answers.txt',
            ),
        )
        write_report(
            os.path.join(
                abs_folder,
                'reports',
                self._tool_name,
                '{0}.txt'.format(self._tool_name),
            ),
            true_answers,
            tool_answers,
            default_ranks,
        )

    def run(self, folder, specification=None, default_ranks=None):
        if default_ranks is None:
            default_ranks = []
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
        self._make_conclusion(abs_folder, default_ranks)

    def name(self):
        return self._tool_name

    def location(self):
        return self._location