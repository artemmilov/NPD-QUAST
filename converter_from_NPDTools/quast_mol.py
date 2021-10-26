import os

from rdkit import Chem


def _is_valid_smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    try:
        Chem.AddHs(m)
    except Exception:
        return False
    return True


class QuastMol:
    def __init__(self, data, scan):
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
        if (self.filename is None) or (self.name is None) or (self.mass is None) or \
                (num_acids is None) or (self.charge is None) or (self.ms_level is None):
            raise Exception
        if not _is_valid_smiles(self.smiles):
            raise Exception
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
