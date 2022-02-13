import os
import rdkit
from subprocess import run, check_output, CalledProcessError

from abstract.abstract_tool import AbstractTool
from quast_mol import QuastMol, QuastMolInitException
from npdtools_database import NpdToolsDatabase


class DereplicatorTool(AbstractTool):
    _spectra_format = 'mgf'
    _database_format = 'csv'
    _tool_name = 'Dereplicator_plus'

    def _init_tool(self, abs_folder):
        super()._init_tool(abs_folder)
        undeployed_database_file = os.path.join(
            abs_folder,
            'temp',
            'database.{0}'.format(
                self._database_format,
            ),
        )
        database_folder = os.path.join(
            abs_folder,
            'temp',
            'tool',
            'deployed_database',
        )
        database = NpdToolsDatabase(database_folder)
        with open(undeployed_database_file, 'r', encoding='utf-8') as undeployed_database:
            mols_data = undeployed_database.read().split('\n\n')
            scan = 1
            for mol_data in mols_data:
                try:
                    quast_mol = QuastMol(mol_data, scan)
                except QuastMolInitException as e:
                    print(e)
                    continue
                database.add_mol(quast_mol)
                scan += 1
        with open(
            os.path.join(
                abs_folder,
                'temp',
                'tool',
                'deployed_database',
                'smiles.info',
            )
        ) as smiles, open(
            os.path.join(
                abs_folder,
                'temp',
                'tool',
                'inches.info',
            ),
            'w',
        ) as inches:
            for line in smiles.readlines():
                try:
                    m = rdkit.Chem.MolFromSmiles(line)
                    inches.write(rdkit.Chem.MolToInchiKey(m))
                except Exception:
                    inches.write('ERROR')

    def _run_tool(self, abs_folder, specification=None):
        path_to_spectres = os.path.join(abs_folder, 'temp', 'spectres')
        path_to_database = os.path.join(abs_folder, 'temp', 'database')
        path_to_results = os.path.join(abs_folder, 'temp', 'results')
        run(
            [
                '/home/artem/Programming/bioinformatics/molDiscovery-2.6.0-beta-Linux/bin/dereplicator.py',
                path_to_spectres,
                '--db-path',
                path_to_database,
                '-o',
                path_to_results,
                '--pass-to-dereplicate',
                '--num_hits_to_report 10',
            ],
        )

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
                                'all_matches.tsv',
                            ),
                    ) as output:
                        for line in output.readlines():
                            answer_id = line.split(' ')[-2][1:-1]
                            answer_inchi_key = list(
                                cur.execute(
                                    'SELECT * FROM molecules WHERE id = {0}' \
                                        .format(answer_id)
                                )
                            )[0][5]
                            tool_answers.write(
                                '{0}${1}\t{2}\t{3}\n'.format(
                                    challenge,
                                    result.split('.')[0],
                                    answer_inchi_key,
                                    str(round(float(line.split(' ')[0]), 3)),
                                ),
                            )
