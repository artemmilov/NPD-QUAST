import os
import shutil
import subprocess

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

from npd_quast.general import NPDQuastError, parse_from_mgf
from npd_quast.tools.abstract_tool import AbstractTool


class _QuastMol:
    def __init__(self, filename, name, mass, smiles):
        self._short_init(filename, name, mass, smiles)

    def _short_init(self, filename, name, mass, smiles):
        self.filename = filename
        self.name = name
        self.mass = mass
        self.smiles = smiles
        if Chem.MolFromSmiles(self.smiles) is None:
            raise NPDQuastError(
                'Smiles {0} is invalid.'.format(self.smiles),
            )
        self.mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))
        self.inchi_key = Chem.inchi.MolToInchiKey(self.mol).split('-')[0]

    def to_mol_file(self, folder):
        Chem.MolToMolFile(self.mol, os.path.join(folder, self.filename), forceV3000=True)


def _check_database_dir(folder):
    if len(os.listdir(folder)) != 3:
        return False
    not_mols = 'mols' not in os.listdir(folder)
    not_library = 'library.info' not in os.listdir(folder)
    not_smiles = 'smiles.info' not in os.listdir(folder)
    if not_mols or not_library or not_smiles:
        return False
    for mol_file in os.listdir(os.path.join(folder, 'mols')):
        if len(mol_file) < 4:
            return False
        if mol_file[-4:] != '.mol':
            return False
    return True


def _rename_database(folder):
    os.rename(
        os.path.join(folder, 'mols'),
        os.path.join(folder, 'old_mols'),
    )
    os.mkdir(os.path.join(folder, 'mols'))
    os.rename(
        os.path.join(folder, 'library.info'),
        os.path.join(folder, 'old_library.info'),
    )
    os.mknod(os.path.join(folder, 'library.info'))
    os.rename(
        os.path.join(folder, 'smiles.info'),
        os.path.join(folder, 'old_smiles.info'),
    )
    os.mknod(os.path.join(folder, 'smiles.info'))


def _clean_database(folder):
    mols = set(os.listdir(os.path.join(folder, 'old_mols')))
    for mol in mols:
        os.remove(os.path.join(folder, 'old_mols', mol))
    os.rmdir(os.path.join(folder, 'old_mols'))
    os.remove(os.path.join(folder, 'old_library.info'))
    os.remove(os.path.join(folder, 'old_smiles.info'))


class _NpdToolsDatabase:
    def __init__(self, folder):
        if not os.path.isdir(folder):
            raise NPDQuastError(
                'There`s no directory {0}.'.format(folder),
            )
        if len(os.listdir(folder)) > 0:
            raise NPDQuastError(
                'Directory {0} is not empty.'.format(folder),
            )
        self.folder = folder
        self._existing_inches = set()
        os.mkdir(os.path.join(folder, 'mols'))
        os.mknod(os.path.join(folder, 'library.info'))
        os.mknod(os.path.join(folder, 'smiles.info'))
        self._l = 0

    def add_mol(self, quast_mol):
        if quast_mol.inchi_key in self._existing_inches:
            return
        quast_mol.to_mol_file(os.path.join(self.folder, 'mols'))
        with open(
                os.path.join(self.folder, 'library.info'),
                'a',
        ) as library_info:
            library_info.write(
                '{0} {1} {2} 1000 DB\n'.format(
                    os.path.join('mols', quast_mol.filename),
                    quast_mol.name,
                    Descriptors.ExactMolWt(quast_mol.mol),
                ),
            )
        with open(
                os.path.join(self.folder, 'smiles.info'),
                'a',
        ) as library_smiles:
            library_smiles.write(Chem.MolToSmiles(quast_mol.mol) + '\n')
        self._l += 1


class AbstractNpdTool(AbstractTool):
    _spectra_format = 'mgf'
    _database_format = 'csv'
    _id_to_inchi = dict()

    def _convert_database(self, from_database, to_database):
        shutil.copyfile(from_database, to_database)

    def _convert_spectra(self, from_spectra, to_spectra):
        shutil.copyfile(from_spectra, to_spectra)

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
        database = _NpdToolsDatabase(database_folder)
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
                    quast_mol = _QuastMol(filename, name, mass, smiles)
                except NPDQuastError:
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
                    m = rdkit.Chem.MolFromSmiles(line)
                    if m is not None:
                        self._id_to_inchi[i] = rdkit.Chem.MolToInchiKey(m).split('-')[0]
                    else:
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
                            str(-float(score)),
                        ),
                    )


class DereplicatorTool(AbstractNpdTool):
    _tool_name = 'Dereplicator'

    def _run_tool(self, abs_folder, specification=None):
        path_to_spectres, path_to_database, path_to_result =\
            super()._run_abstract_tool(abs_folder, specification)
        subprocess.run(
            [
                self._location,
                path_to_spectres,
                '--db-path',
                path_to_database,
                '-o',
                path_to_result,
            ],
            capture_output=True,
        )


class DereplicatorPlusTool(AbstractNpdTool):
    _tool_name = 'Dereplicator_plus'

    def _run_tool(self, abs_folder, specification=None):
        path_to_spectres, path_to_database, path_to_result =\
            super()._run_abstract_tool(abs_folder, specification)
        subprocess.run(
            [
                self._location,
                path_to_spectres,
                '--db-path',
                path_to_database,
                '-o',
                path_to_result,
                '--pass-to-dereplicate',
                '--num_hits_to_report 100',
                '--min-score',
                '1',
            ],
            capture_output=True,
        )
