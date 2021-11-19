import os

from rdkit import Chem


def _is_valid_smiles(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
    except TypeError:
        return False
    if m is None:
        return False
    return True


class QuastMolInitException(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return 'QuastMolInitException: {0}'.format(self.message)


class QuastMol:
    def __init__(self, *args):
        if len(args) == 2:
            data, scan = args
            self._full_init(data, scan)
        elif len(args) == 4:
            filename, name, mass, smiles = args
            self._short_init(filename, name, mass, smiles)
        else:
            raise QuastMolInitException(
                'Expected 2 or 4 arguments, but {0} given.'.format(
                    len(args),
                )
            )

    def _full_init(self, data, scan):
        self.scan = scan
        self.spectra = ''
        self.smiles, self.filename, self.name, self.mass, num_acids = None, None, None, None, None
        self.charge, self.ms_level = None, None
        self.filename = f'compound_{self.scan:05d}.mol'
        num_acids = 1000
        for line in data.splitlines():
            if '=' in line:
                line_split = line.split('=', maxsplit=1)
                left, right = line_split[0], line_split[1]
                if (left == 'SMILES') and ('.' not in right):
                    self.smiles = right
                elif (left == 'NAME') and (left != ''):
                    self.name = right.replace(' ', '_')
                elif left == 'PEPMASS':
                    self.mass = float(right)
                elif left == 'CHARGE':
                    self.charge = right
                elif left == 'MSLEVEL':
                    self.ms_level = right
            elif (line != 'BEGIN IONS') and (line != 'END IONS') and (line != ''):
                self.spectra += line + '\n'
        unfounded_params = []
        if self.filename is None:
            unfounded_params.append('filename')
        if self.name is None:
            unfounded_params.append('name')
        if self.mass is None:
            unfounded_params.append('mass')
        if num_acids is None:
            unfounded_params.append('num_acids')
        if self.charge is None:
            unfounded_params.append('charge')
        if self.ms_level is None:
            unfounded_params.append('ms_level')
        if len(unfounded_params) == 1:
            raise QuastMolInitException('Scan: {0}. {1} is unfounded.'.format(scan, unfounded_params[0]))
        elif len(unfounded_params) == 2:
            raise QuastMolInitException('Scan: {0}. {1} are unfounded.'.format(scan, ','.join(unfounded_params)))
        if not _is_valid_smiles(self.smiles):
            raise QuastMolInitException('Scan: {0}. Smiles {1} is invalid.'.format(scan, self.smiles))
        self.mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))
        self.inchi_key = Chem.inchi.MolToInchiKey(self.mol).split('-')[0]

    def _short_init(self, filename, name, mass, smiles):
        self.filename = filename
        self.name = name
        self.mass = mass
        self.smiles = smiles
        if not _is_valid_smiles(self.smiles):
            raise QuastMolInitException('Smiles {0} is invalid.'.format(self.smiles))
        self.mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))
        self.inchi_key = Chem.inchi.MolToInchiKey(self.mol).split('-')[0]

    def to_npdtools_input(self):
        return '''BEGIN IONS\nSCANS={0}\nPEPMASS={1}\nCHARGE={2}\nMSLEVEL={3}\n{4}END IONS\n'''.format(
            self.scan,
            self.mass,
            self.charge,
            self.ms_level,
            self.spectra,
        )

    def to_mol_file(self, folder):
        Chem.MolToMolFile(self.mol, os.path.join(folder, self.filename), forceV3000=True)
