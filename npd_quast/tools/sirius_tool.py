import os
import shutil
import subprocess

from npd_quast.general import parse_from_mgf
from npd_quast.tools.abstract_tool import AbstractTool


class SiriusTool(AbstractTool):
    _spectra_format = 'mgf'
    _database_format = 'txt'
    _tool_name = 'Sirius'

    def _init_tool(self, abs_folder):
        super()._init_tool(abs_folder)
        os.mkdir(
            os.path.join(abs_folder, 'temp', 'tool', 'cur_results'),
        )

    def _convert_database(self, from_database, to_database):
        with open(from_database) as csv, \
                open(to_database, 'w') as txt:
            for line in csv.readlines()[1:-1]:
                identifier, name, _, _, smiles, _, _ = parse_from_mgf(line)
                txt.write('{0}\t{1}\t{2}\n'.format(smiles, identifier, name))

    def _convert_spectra(self, from_spectra, to_spectra):
        shutil.copyfile(from_spectra, to_spectra)

    def _run_tool(self, abs_folder, specification=None):
        path_to_spectres = os.path.join(abs_folder, 'temp', 'spectres')
        path_to_results = os.path.join(abs_folder, 'temp', 'tool', 'cur_results')
        path_to_database = os.path.join(abs_folder, 'temp', 'database.txt')
        if os.path.isdir(path_to_results):
            shutil.rmtree(path_to_results)
            os.mkdir(path_to_results)

        for spectra in os.listdir(path_to_spectres):
            try:
                subprocess.run(
                    '''export PATH=\"{0}:$PATH\"; \
sirius -i {1} custom-db --name cur_database; \
sirius -i {2} -o {3} formula -c 10 structure --database cur_database'''.format(
                        self._location,
                        path_to_database,
                        os.path.join(path_to_spectres, spectra),
                        os.path.join(path_to_results, spectra.split('.')[0]),
                    ),
                    shell=True,
                    capture_output=True,
                )
            except subprocess.CalledProcessError:
                pass

    def _parse_output(self, abs_folder, challenge_name):
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
                        self._tool_name,
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
                            '{0}${1}\t{2}\t{3}\n'.format(
                                challenge_name,
                                result,
                                answer_inchi_key,
                                str(-float(score)),
                            ),
                        )
