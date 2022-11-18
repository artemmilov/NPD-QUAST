import os
import shutil
import subprocess

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors


def parse_from_mgf(s):
    res = []
    cur = ''
    in_brackets = False
    for a in s:
        if (a == ',') and (not in_brackets):
            res.append(cur)
            cur = ''
        elif a == '\"':
            in_brackets = not in_brackets
        else:
            cur += a
    if cur != '':
        res.append(cur)
    return res


class _QuastMol:
    def __init__(self, filename, name, mass, smiles):
        self._short_init(filename, name, mass, smiles)

    def _short_init(self, filename, name, mass, smiles):
        self.filename = filename
        self.name = name
        self.mass = mass
        self.smiles = smiles
        if Chem.MolFromSmiles(self.smiles) is None:
            raise Exception(
                'Smiles {0} is invalid.'.format(self.smiles),
            )
        self.mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))
        self.inchi_key = Chem.inchi.MolToInchiKey(self.mol).split('-')[0]

    def to_mol_file(self, folder):
        Chem.MolToMolFile(self.mol, os.path.join(folder, self.filename), forceV3000=True)


class _NpdToolsDatabase:
    def __init__(self):
        self._existing_inches = set()
        os.mkdir(snakemake.output[0])
        os.mkdir(os.path.join(snakemake.output[0], 'mols'))
        open(os.path.join(snakemake.output[0], 'library.info'), 'w')
        open(os.path.join(snakemake.output[0], 'library.smiles'), 'w')
        self._l = 0

    def add_mol(self, quast_mol):
        if quast_mol.inchi_key in self._existing_inches:
            return
        quast_mol.to_mol_file(os.path.join(snakemake.output[0], 'mols'))
        with open(
                os.path.join(snakemake.output[0], 'library.info'),
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
                os.path.join(snakemake.output[0], 'library.smiles'),
                'a',
        ) as library_smiles:
            library_smiles.write(Chem.MolToSmiles(quast_mol.mol) + '\n')
        self._l += 1


def deploy_database():
    # with open("challenges/Challenge-082/database.csv", 'r', encoding='utf-8') as raw_database:
    with open(snakemake.input[0], 'r', encoding='utf-8') as raw_database:
        database = _NpdToolsDatabase()
        mols_data = raw_database.readlines()
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
            except Exception:
                continue
            database.add_mol(quast_mol)
            scan += 1


deploy_database()
# open("temp/dereplicator+/{0}/".format(os.path.split(snakemake.output[0])[1]), 'w')
# open("temp/dereplicator+/Challenge-082/", 'w')
