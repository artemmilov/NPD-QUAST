import os

from rdkit import Chem


class AbstractMol:
    filename = None
    mol = None

    def _full_init(self, *args):
        pass

    def _short_init(self, *args):
        pass

    def to_tool_input(self, *args):
        pass

    def to_mol_file(self, folder):
        Chem.MolToMolFile(self.mol, os.path.join(folder, self.filename), forceV3000=True)
