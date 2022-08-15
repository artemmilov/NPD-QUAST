import configparser
import os
import shutil


class AbstractTool:
    _specter_format = None
    _database_format = None
    _tool_name = None
    _location = None

    def __init__(self):
        config = configparser.ConfigParser()
        config.read('npd_quast.ini')
        self._location = self.clarify_location(config['supported tools'][self._tool_name])

    def _init_tool(self, abs_folder, report):
        if not os.path.isdir(os.path.join(abs_folder, 'temp')):
            os.mkdir(os.path.join(abs_folder, 'temp'))
        if len(os.listdir(os.path.join(abs_folder, 'temp'))) > 0:
            for subdir in os.listdir(
                    os.path.join(abs_folder, 'temp'),
            ):
                shutil.rmtree(subdir)
        os.mkdir(os.path.join(abs_folder, 'temp', 'spectra'))
        os.mkdir(os.path.join(abs_folder, 'temp', 'tool'))
        if report in os.listdir(
                os.path.join(abs_folder, 'reports'),
        ):
            shutil.rmtree(
                os.path.join(abs_folder, 'reports', report),
            )
        os.mkdir(
            os.path.join(abs_folder, 'reports', report),
        )
        os.mknod(
            os.path.join(
                abs_folder,
                'reports',
                report,
                'tool_answers.txt',
            ),
        )

    def clarify_location(self, loc):
        pass

    def _convert_database(self, from_database, to_database):
        pass

    def _write_database(self, abs_folder, database):
        if os.path.exists(os.path.join(
                abs_folder,
                'temp',
                'database.{0}'.format(self._database_format),
        )):
            os.remove(os.path.join(
                abs_folder,
                'temp',
                'database.{0}'.format(self._database_format),
            ))
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

    def _convert_specter(self, from_specter, to_specter):
        pass

    def _write_spectra(self, abs_folder, spectra):
        if os.path.isdir(os.path.join(abs_folder, 'temp', 'spectra')):
            shutil.rmtree(os.path.join(abs_folder, 'temp', 'spectra'))
            os.mkdir(os.path.join(abs_folder, 'temp', 'spectra'))
        if not os.path.isdir(spectra):
            raise ValueError('{0} is not a directory'.format(spectra))
        for specter in os.listdir(spectra):
            if len(specter.split('.')) != 2:
                raise ValueError('Indefinable format')
            self._convert_specter(
                os.path.join(spectra, specter),
                os.path.join(abs_folder, 'temp', 'spectra', specter),
            )

    def _run_tool(self, abs_folder, specification, logger):
        pass

    def _parse_output(self, abs_folder, challenge_name, report):
        pass

    def _run_challenge(self, abs_folder, challenge, report, specification, logger, debug):
        if debug:
            challenge_name = os.path.split(challenge)[-1]
            os.mkdir(os.path.join(abs_folder, 'debug', report, challenge_name))
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
        self._write_spectra(
            abs_folder,
            os.path.join(challenge, 'spectra'),
        )
        completed_process = self._run_tool(abs_folder, specification, logger)
        self._parse_output(abs_folder, os.path.split(challenge)[-1], report)
        if debug:
            challenge_name = os.path.split(challenge)[-1]
            shutil.copytree(
                os.path.join(abs_folder, 'temp'),
                os.path.join(abs_folder, 'debug', report, challenge_name, 'work')
            )
            with open(
                    os.path.join(abs_folder, 'debug', report, challenge_name, 'stdout'),
                    'w'
            ) as stdout:
                stdout.write(completed_process.stdout.decode())
            with open(
                    os.path.join(abs_folder, 'debug', report, challenge_name, 'stderr'),
                    'w'
            ) as stderr:
                stderr.write(completed_process.stderr.decode())

    def run(self, folder, report, specification, logger, debug=False):
        abs_folder = os.path.abspath(folder)
        logger.info('Initialization...')
        self._init_tool(abs_folder, report)
        for challenge in os.listdir(
                os.path.join(abs_folder, 'challenges'),
        ):
            logger.info('Run challenge {}...'.format(challenge))
            self._run_challenge(
                abs_folder,
                os.path.join(abs_folder, 'challenges', challenge),
                report,
                specification,
                logger,
                debug,
            )
            logger.info('Ok!')

    def name(self):
        return self._tool_name

    def location(self):
        return self._location
