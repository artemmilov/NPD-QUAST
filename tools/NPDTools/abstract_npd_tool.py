import os
import shutil

import rdkit
from tools.NPDTools.npdtools_database import NpdToolsDatabase
from tools.NPDTools.quast_mol import QuastMol, QuastMolInitException

from tools.abstract_tool import AbstractTool
from general import parse_from_mgf


class AbstractNpdTool(AbstractTool):
    _spectra_format = 'mgf'
    _database_format = 'csv'
    _id_to_inchi = dict()

    def _deploy_database(self, abs_folder):
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
        if os.path.isdir(database_folder):
            shutil.rmtree(database_folder)
        os.mkdir(database_folder)
        database = NpdToolsDatabase(database_folder)
        with open(undeployed_database_file, 'r', encoding='utf-8') as undeployed_database:
            mols_data = undeployed_database.readlines()
            scan = 1
            for mol_data in mols_data:
                try:
                    filename = os.path.join(
                        parse_from_mgf(mol_data)[0] + '.mol',
                    )
                    name = parse_from_mgf(mol_data)[0]
                    mass = parse_from_mgf(mol_data)[2]
                    smiles = parse_from_mgf(mol_data)[4]
                    quast_mol = QuastMol(filename, name, mass, smiles)
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
        ) as smiles:
            for i, line in enumerate(smiles.readlines()):
                if line != '':
                    try:
                        m = rdkit.Chem.MolFromSmiles(line)
                        self._id_to_inchi[i] = rdkit.Chem.MolToInchiKey(m).split('-')[0]
                    except Exception:
                        self._id_to_inchi[i] = 'ERROR'

    def _run_abstract_tool(self, abs_folder, specification=None):
        self._deploy_database(abs_folder)
        path_to_spectres = os.path.join(abs_folder, 'temp', 'spectres')
        path_to_database = os.path.join(abs_folder, 'temp', 'tool', 'deployed_database')
        path_to_result = os.path.join(abs_folder, 'temp', 'tool', 'cur_result')
        if os.path.isdir(path_to_result):
            shutil.rmtree(path_to_result)
        os.mkdir(path_to_result)
        return path_to_spectres, path_to_database, path_to_result

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
            with open(
                    os.path.join(
                        abs_folder,
                        'temp',
                        'tool',
                        'cur_result',
                        'all_matches.tsv',
                    ),
            ) as output:
                for line in output.readlines()[1:]:
                    answer_id = int(line.split('\t')[3])
                    answer_inchi_key = self._id_to_inchi[answer_id]
                    spectra = os.path.split(line.split('\t')[0])[-1].split('.')[0]
                    score = line.split('\t')[5]
                    tool_answers.write(
                        '{0}${1}\t{2}\t{3}\n'.format(
                            challenge_name,
                            spectra,
                            answer_inchi_key,
                            score,
                        ),
                    )
