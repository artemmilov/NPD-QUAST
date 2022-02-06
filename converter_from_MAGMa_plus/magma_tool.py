import os
from subprocess import check_output

from abstract.abstract_tool import AbstractTool


class MagmaTool(AbstractTool):
    _spectra_format = 'tree'
    _database_format = 'db'
    _tool_name = 'MAGMa_plus'

    def _init_tool(self, abs_folder):
        super()._init_tool(abs_folder)
        with open(
                os.path.join(
                    abs_folder,
                    'temp',
                    'tool',
                    'magma_job.ini',
                ),
                'w',
        ) as magma_ini:
            magma_ini.write(
                '''[magma job]

# Location of structure database to fetch candidate molecules to match against ms peak trees
structure_database.hmdb = {0}

chemical_engine = rdkit'''.format(
                    os.path.join(
                        abs_folder,
                        'temp',
                        'database.db',
                    ),
                ),
            )
        with open(
                os.path.join(
                    abs_folder,
                    'temp',
                    'tool',
                    'script.txt',
                ),
                'w',
        ) as script:
            script.write('''#!/bin/bash
cd {0}
export PATH=/home/artem/Programming/miniconda3/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate magma-plus-env
export MAGMAPLUS_CLASSIFIER_PATH=/home/artem/Programming/bioinformatics/MAGMa-plus
path_to_magma=/home/artem/Programming/bioinformatics/MAGMa-plus/MAGMa_plus.py
python $path_to_magma read_ms_data -i 1 -p 5 -q 0.001 -f mass_tree $1 $2
python $path_to_magma annotate -c 0 -d 0 -b 3 -w 1 -s hmdb $2
python $path_to_magma export_result $2'''.format(
                os.path.join(abs_folder, 'temp', 'tool'),
            ),
            )
        with open(
                os.path.join(
                    abs_folder,
                    'temp',
                    'tool',
                    'structure.db',
                ),
                'w',
        ) as script:
            script.write(
                '{0}'.format(
                    os.path.join(
                        abs_folder,
                        'temp',
                        'database.db',
                    ),
                ),
            )

    def _run_tool(self, abs_folder, specification=None):
        path_to_script = os.path.join(abs_folder, 'temp', 'tool', 'script.txt')
        path_to_spectres = os.path.join(abs_folder, 'temp', 'spectres')
        path_to_results = os.path.join(abs_folder, 'temp', 'results')
        for spectra in os.listdir(path_to_spectres):
            result = spectra.split('.')[0] + '.db'
            output = check_output(
                ' '.join(
                    [
                        'bash',
                        path_to_script,
                        os.path.join(path_to_spectres, spectra),
                        os.path.join(path_to_results, result.split('.')[0]),
                    ],
                ),
                shell=True,
            ).decode('utf-8')
            os.mkdir(
                os.path.join(
                    abs_folder, 'temp', 'results', result.split('.')[0], 'data',
                )
            )
            with open(
                    os.path.join(
                        abs_folder, 'temp', 'results', result.split('.')[0], 'data', result.split('.')[0] + '.txt',
                    ),
                    'w',
            ) as res:
                write = False
                for line in output.splitlines():
                    if (line == '') or (not line[0].isnumeric()):
                        write = False
                    if write:
                        res.write(line + '\n')
                    if line == 'Candidate_Score Name Smiles':
                        write = True

    def _parse_output(self, abs_folder):
        if not os.path.isdir(
            os.path.join(
                abs_folder,
                'reports',
                self._tool_name,
            ),
        ):
            os.mkdir(
                os.path.join(
                    abs_folder,
                    'reports',
                    self._tool_name,
                ),
            )

        with open(
                os.path.join(
                    abs_folder,
                    'reports',
                    self._tool_name,
                    'tool_answers.txt',
                ),
                'w',
        ) as tool_answers:
            for challenge in filter(
                    lambda c: os.path.isdir(
                        os.path.join(
                            abs_folder,
                            'temp',
                            'results',
                            c
                        )
                    ),
                    os.listdir(
                        os.path.join(abs_folder, 'temp', 'results'),
                    ),
            ):
                for result in os.listdir(
                        os.path.join(
                            abs_folder,
                            'temp',
                            'results',
                            challenge,
                            'data'
                        ),
                ):
                    with open(
                            os.path.join(
                                abs_folder,
                                'temp',
                                'results',
                                challenge,
                                'data',
                                result,
                            ),
                    ) as output:
                        for line in output.readlines():
                            tool_answers.write(
                                '{0}${1}\t{2}\t{3}\n'.format(
                                    challenge,
                                    result.split('.')[0],
                                    line.split(' ')[-2][1:-1],
                                    line.split(' ')[0],
                                ),
                            )
