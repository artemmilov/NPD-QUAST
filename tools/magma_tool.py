import os
import shutil
import sqlite3
import zlib
from subprocess import CalledProcessError, check_output

from rdkit import Chem
from rdkit.Chem import Crippen

from general import parse_from_mgf, NPDQuastError
from tools.abstract_tool import AbstractTool


def _parse_line_for_magma(line):
    parsed_line = parse_from_mgf(line)
    reference = scan = parsed_line[0]
    name = parsed_line[1]
    mass = float(parsed_line[2])
    mol_form = parsed_line[3]
    smiles = parsed_line[4]
    _ = parsed_line[5]
    inchi_key = parsed_line[6]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise NPDQuastError('Wrong smiles: {0}'.format(smiles))
    return (
        scan,
        int(round(mass * 1e6)),
        Chem.GetFormalCharge(mol),
        mol.GetNumHeavyAtoms(),
        sqlite3.Binary(
            zlib.compress(
                Chem.MolToMolBlock(mol).encode('utf-8'),
            ),
        ),
        inchi_key.split('-')[0],
        mol_form,
        name,
        reference,
        int(round(Crippen.MolLogP(mol) * 10)),
    )


def _create_sqlite_database(
        raw_database_file,
        converted_database_file,
):
    conn = sqlite3.connect(converted_database_file)
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE molecules (id TEXT PRIMARY KEY,
                                             mim INTEGER NOT NULL,
                                             charge INTEGER NOT NULL,
                                             natoms INTEGER NOT NULL,
                                             molblock TEXT,
                                             inchikey TEXT,
                                             molform TEXT,
                                             name TEXT,
                                             reference TEXT,
                                             logp INT)''')
        conn.commit()
    except:
        return

    with open(raw_database_file) as raw_database:
        for line in raw_database.readlines()[1:]:
            if line != '':
                try:
                    c.execute(
                        '''INSERT INTO molecules (id, mim, charge, natoms, molblock, inchikey,
                         molform, name, reference, logp) VALUES (?,?,?,?,?,?,?,?,?,?)''',
                        _parse_line_for_magma(line),
                    )
                except NPDQuastError:
                    pass
    conn.commit()

    c.execute('PRAGMA temp_store = 2')
    c.execute('''CREATE INDEX idx_cover ON molecules (charge, mim, natoms, reference, 
molform, inchikey, name, molblock, logp)''')
    conn.commit()


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
export PATH={1}:$PATH
eval "$(conda shell.bash hook)"
conda activate magma-plus-env
export MAGMAPLUS_CLASSIFIER_PATH=/home/artem/Programming/bioinformatics/MAGMa-plus
path_to_magma=/home/artem/Programming/bioinformatics/MAGMa-plus/MAGMa_plus.py
python $path_to_magma read_ms_data -i 1 -p 5 -q 0.001 -f mass_tree $1 $2
python $path_to_magma annotate -c 0 -d 0 -b 3 -w 1 -s hmdb $2
python $path_to_magma export_result $2'''.format(
                    os.path.join(abs_folder, 'temp', 'tool'),
                    self._location,
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

    def _convert_database(self, from_database, to_database):
        _create_sqlite_database(from_database, to_database)

    def _convert_spectra(self, from_spectra, to_spectra):
        with open(from_spectra) as mgf:
            mzs = []
            norm_intensities = []
            mass = None
            for line in mgf.readlines():
                if '=' in line:
                    if line.split('=')[0] == 'PEPMASS':
                        mass = float(line.split('=')[1])
                elif len(line.split('\t')) == 2:
                    mzs.append(float(line.split('\t')[0]))
                    norm_intensities.append(float(line.split('\t')[1]))
            normalizer = max(norm_intensities)
            record = '{0}: 100 ('.format(mass)
            record += ', '.join(
                '{0}: {1}'.format(
                    mzs[i],
                    norm_intensities[i] * 100 / normalizer,
                )
                for i in range(0, len(mzs))
            )
            record += ')'
            with open(to_spectra, 'w') as mass_tree:
                mass_tree.write(record)

    def _run_tool(self, abs_folder, specification=None):
        my_env = os.environ.copy()
        my_env["PATH"] = "/home/artem/Programming/miniconda3/bin:" + my_env["PATH"]
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
                    env=my_env,
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
