import os
import shutil
import sqlite3
from subprocess import check_output, CalledProcessError

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
        os.mkdir(
            os.path.join(
                abs_folder,
                'temp',
                'tool',
                'cur_results',
            ),
        )
        os.mkdir(
            os.path.join(
                abs_folder,
                'temp',
                'tool',
                'cur_trees',
            ),
        )

    def _run_tool(self, abs_folder, specification=None):
        path_to_script = os.path.join(abs_folder, 'temp', 'tool', 'script.txt')
        path_to_spectres = os.path.join(abs_folder, 'temp', 'spectres')
        path_to_trees = os.path.join(abs_folder, 'temp', 'tool', 'cur_trees')
        path_to_results = os.path.join(abs_folder, 'temp', 'tool', 'cur_results')
        if os.path.isdir(path_to_trees):
            shutil.rmtree(path_to_trees)
            os.mkdir(path_to_trees)
        if os.path.isdir(path_to_results):
            shutil.rmtree(path_to_results)
            os.mkdir(path_to_results)
        for spectra in os.listdir(path_to_spectres):
            converted_tree = spectra.split('.')[0] + '.db'
            try:
                output = check_output(
                    [
                        'bash',
                        path_to_script,
                        os.path.join(path_to_spectres, spectra),
                        os.path.join(path_to_trees, converted_tree),
                    ],
                ).decode('utf-8')
            except CalledProcessError:
                output = ''
            with open(
                    os.path.join(
                        path_to_results,
                        converted_tree.split('.')[0] + '.txt',
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

    def _parse_output(self, abs_folder, challenge_name):
        with open(
                os.path.join(
                    abs_folder,
                    'reports',
                    self._tool_name,
                    'tool_answers.txt',
                ),
                'a',
        ) as tool_answers:
            for result in os.listdir(
                    os.path.join(abs_folder, 'temp', 'tool', 'cur_results'),
            ):
                conn = sqlite3.connect(
                    os.path.join(
                        abs_folder,
                        'temp',
                        'database.db',
                    )
                )
                cur = conn.cursor()
                with open(
                        os.path.join(
                            abs_folder,
                            'temp',
                            'tool',
                            'cur_results',
                            result,
                        ),
                ) as output:
                    for line in output.readlines():
                        answer_id = line.split(' ')[-2][1:-1]
                        answer_inchi_key = list(
                            cur.execute(
                                'SELECT * FROM molecules WHERE id = {0}'.format(
                                    answer_id,
                                ),
                            )
                        )[0][5]
                        tool_answers.write(
                            '{0}${1}\t{2}\t{3}\n'.format(
                                challenge_name,
                                result.split('.')[0],
                                answer_inchi_key,
                                str(round(float(line.split(' ')[0]), 3)),
                            ),
                        )
