import os
import shutil
import subprocess

from npd_quast.general import parse_from_mgf
from npd_quast.tools.abstract_tool import AbstractTool


class SiriusTool(AbstractTool):
    _specter_format = 'mgf'
    _database_format = 'txt'
    _tool_name = 'Sirius'

    def _init_tool(self, abs_folder, report):
        super()._init_tool(abs_folder, report)
        os.mkdir(
            os.path.join(abs_folder, 'temp', 'tool', 'cur_results'),
        )

    def _convert_database(self, from_database, to_database):
        with open(from_database) as csv, \
                open(to_database, 'w') as txt:
            for line in csv.readlines()[1:-1]:
                identifier, name, _, _, smiles, _, _ = parse_from_mgf(line)
                txt.write('{0}\t{1}\t{2}\n'.format(smiles, identifier, name))

    def _convert_specter(self, from_specter, to_specter):
        shutil.copyfile(from_specter, to_specter)

    def _run_tool(self, abs_folder, specification, logger):
        path_to_spectra = os.path.join(abs_folder, 'temp', 'spectra')
        path_to_results = os.path.join(abs_folder, 'temp', 'tool', 'cur_results')
        path_to_database = os.path.join(abs_folder, 'temp', 'database.txt')
        if os.path.isdir(path_to_results):
            shutil.rmtree(path_to_results)
            os.mkdir(path_to_results)
        for specter in os.listdir(path_to_spectra):
            command = 'export PATH=\"{0}:$PATH\";'.format(self._location)
            command += 'sirius -i {0} {1} custom-db {2};'.format(
                path_to_database,
                ' '.join(
                    [
                        '{0} {1}'.format(k, v)
                        for k, v in specification.items()
                        if (v is not None) and (not isinstance(v, dict))
                    ] +
                    [
                        '{0}'.format(k, v)
                        for k, v in specification.items()
                        if (v is None) and (not isinstance(v, dict))
                    ]
                ),
                ' '.join(
                    [
                        '{0} {1}'.format(k, v)
                        for k, v in specification['custom-db'].items()
                        if (v is not None) and (not isinstance(v, dict))
                    ] +
                    [
                        '{0}'.format(k, v)
                        for k, v in specification['custom-db'].items()
                        if (v is None) and (not isinstance(v, dict))
                    ]
                ),
            )
            command += 'sirius -i {0} -o {1} {2} formula {3} structure {4}'.format(
                os.path.join(path_to_spectra, specter),
                os.path.join(path_to_results, specter.split('.')[0]),
                ' '.join(
                    [
                        '{0} {1}'.format(k, v)
                        for k, v in specification.items()
                        if (v is not None) and (not isinstance(v, dict))
                    ] +
                    [
                        '{0}'.format(k, v)
                        for k, v in specification.items()
                        if (v is None) and (not isinstance(v, dict))
                    ]
                ),
                ' '.join(
                    [
                        '{0} {1}'.format(k, v)
                        for k, v in specification["formula"].items()
                        if (v is not None) and (not isinstance(v, dict))
                    ] +
                    [
                        '{0}'.format(k, v)
                        for k, v in specification["formula"].items()
                        if (v is None) and (not isinstance(v, dict))
                    ]
                ),
                ' '.join(
                    [
                        '{0} {1}'.format(k, v)
                        for k, v in specification["formula"]["structure"].items()
                        if (v is not None) and (not isinstance(v, dict))
                    ] +
                    [
                        '{0}'.format(k, v)
                        for k, v in specification["formula"]["structure"].items()
                        if (v is None) and (not isinstance(v, dict))
                    ]
                ),
            )
            logger.info('Run commands:\n{}'.format(command.replace(';', '\n')))
            try:
                completed_process = subprocess.run(
                    command,
                    shell=True,
                    capture_output=True,
                )
            except subprocess.CalledProcessError as e:
                logger.error(e)
                return
            return completed_process

    def _parse_output(self, abs_folder, challenge_name, report):
        for result in os.listdir(
                os.path.join(
                    abs_folder,
                    'temp',
                    'tool',
                    'cur_results',
                ),
        ):
            with open(
                    os.path.join(
                        abs_folder,
                        'reports',
                        report,
                        'tool_answers.txt',
                    ),
                    'a',
            ) as tool_answers:
                with open(
                        os.path.join(
                            abs_folder,
                            'temp',
                            'tool',
                            'cur_results',
                            result,
                            '0_{0}_FEATURE_1'.format(result),
                            'structure_candidates.tsv',
                        ),
                ) as output:
                    for line in output.readlines()[1:]:
                        answer_inchi_key = line.split('\t')[5]
                        score = line.split('\t')[2]
                        tool_answers.write(
                            '{0}\t{1}\t{2}\t{3}\n'.format(
                                challenge_name,
                                result,
                                answer_inchi_key,
                                str(-float(score)),
                            ),
                        )
