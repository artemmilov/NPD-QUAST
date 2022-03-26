import os

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi

from abstract.abstract_database import AbstractDatabase
from last.NPDTools.quast_mol import QuastMol, QuastMolInitException
from general import parse_from_mgf


class DatabaseInitException(Exception):
    pass


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


class NpdToolsDatabase(AbstractDatabase):
    def _empty_init(self, folder):
        self.folder = folder
        self._existing_inches = set()
        os.mkdir(os.path.join(folder, 'mols'))
        os.mknod(os.path.join(folder, 'library.info'))
        os.mknod(os.path.join(folder, 'smiles.info'))
        self._l = 0

    def _not_empty_init(self, folder):
        self.folder = folder
        self._existing_inches = set()
        self._l = 0
        if not _check_database_dir(folder):
            raise DatabaseInitException(
                'Dir {0} is not database'.format(folder),
            )
        _rename_database(folder)
        with open(os.path.join(folder, 'old_library.info')) as library_info, \
                open(os.path.join(folder, 'old_smiles.info')) as smiles_info:
            for line, smiles in list(
                    zip(
                        library_info.readlines(),  # [:-1],
                        smiles_info.readlines(),  # [:-1],
                    ),
            ):
                if len(line.split()) != 5:
                    raise DatabaseInitException()
                try:
                    quast_mol = QuastMol(
                        os.path.split(line.split()[0])[-1],
                        line.split()[1],
                        float(line.split()[2]),
                        smiles,
                    )
                except QuastMolInitException:
                    continue
                self.add_mol(quast_mol)
        _clean_database(folder)

    def clean(self):
        os.remove(
            os.path.join(self.folder, 'library.info')
        )
        os.remove(
            os.path.join(self.folder, 'smiles.info')
        )
        for mol in os.listdir(
                os.path.join(self.folder, 'mols')
        ):
            os.remove(
                os.path.join(self.folder, 'mols', mol)
            )
        os.rmdir(
            os.path.join(self.folder, 'mols')
        )
        self._empty_init(self.folder)

    def transform_from_casmi_data(self, casmi_file):
        folder = os.path.dirname(casmi_file)
        name = os.path.split(casmi_file)[-1].replace(' ', '')
        if len(name.split('.')) != 2:
            raise DatabaseInitException(
                'Incorrect name: {0}'.format(name)
            )
        if name.split('.')[1] != 'csv':
            raise DatabaseInitException()
        if name.split('.')[0] in os.listdir(folder):
            raise DatabaseInitException()
        with open(casmi_file) as casmi_data:
            for line in casmi_data.readlines()[1:]:
                if line != '':
                    splitted_line = parse_from_mgf(line)
                    identifier = int(line.split(',')[0])
                    name = splitted_line[1]
                    mass = float(splitted_line[2])
                    smiles = splitted_line[4]
                    try:
                        self.add_mol(
                            QuastMol(
                                f'compound_{identifier:09d}.mol',
                                name,
                                mass,
                                smiles,
                            ),
                        )
                    except QuastMolInitException:
                        pass

    def __getitem__(self, i):
        with open(os.path.join(self.folder, 'library.info')) as library_info, \
                open(os.path.join(self.folder, 'smiles.info')) as library_smiles:
            mol_data = library_info.readlines()[i]
            full_filename, name, mass = mol_data.split(' ')[:3]
            smiles = library_smiles.readlines()[i]
            return QuastMol(
                os.path.split(full_filename)[-1],
                name,
                mass,
                smiles,
            )

    def add_mol(self, quast_mol):
        if quast_mol.inchi_key in self._existing_inches:
            return False
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
        return True


def _check_database_folder(database_folder):
    if len(os.listdir(database_folder)) != 3:
        return False
    if 'library.info' not in os.listdir(database_folder):
        return False
    if 'smiles.info' not in os.listdir(database_folder):
        return False
    if 'mols' not in os.listdir(database_folder):
        return False
    return True


def merge_databases(database_folder_1, database_folder_2, result_database_folder):
    os.mkdir(os.path.join(result_database_folder, 'mols'))
    with open(os.path.join(result_database_folder, 'library.info'), 'w') as result_library_info, \
            open(os.path.join(result_database_folder, 'smiles.info'), 'w') as result_library_smiles:
        existing_inches = set()
        with open(os.path.join(database_folder_1, 'library.info')) as library_info_1, \
                open(os.path.join(database_folder_1, 'smiles.info')) as library_smiles_1, \
                open(os.path.join(database_folder_2, 'library.info')) as library_info_2, \
                open(os.path.join(database_folder_2, 'smiles.info')) as library_smiles_2:
            for i, database in enumerate(
                    [
                        list(zip(library_info_1.readlines()[:-1], library_smiles_1.readlines()[:-1])),
                        list(zip(library_info_2.readlines()[:-1], library_smiles_2.readlines()[:-1])),
                    ],
            ):
                for line, smiles in database:
                    mol_file = line.split(' ')[0]
                    if i == 0:
                        m = Chem.MolFromMolFile(os.path.join(database_folder_1, mol_file))
                    else:
                        m = Chem.MolFromMolFile(os.path.join(database_folder_2, mol_file))
                    if m is not None:
                        if inchi.MolToInchiKey(m).split('.')[0] not in existing_inches:
                            existing_inches.add(inchi.MolToInchiKey(m).split('.')[0])
                            result_library_info.write(line)
                            result_library_smiles.write(smiles)
                            Chem.MolToMolFile(
                                m,
                                os.path.join(result_database_folder, mol_file),
                                forceV3000=True,
                            )


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    database_folder_1, database_folder_2, result_database_folder = input().split(' ')
    database_1 = NpdToolsDatabase(database_folder_1)
    database_2 = NpdToolsDatabase(database_folder_2)
    result_database = NpdToolsDatabase(result_database_folder)
    result_database.add_database(database_1)
    result_database.add_database(database_2)


if __name__ == '__main__':
    main()
    # ../files/output_sample/database ../files/pnpdatabase ../files/merged_database
