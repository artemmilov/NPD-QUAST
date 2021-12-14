from rdkit import Chem
from rdkit.Chem import Crippen

from abstract.abstract_mol import AbstractMol
from general import parse_from_mgf


class MagmaPlusInitException(Exception):
    pass


class MagmaPlusMol(AbstractMol):
    def __init__(self, *args):
        if len(args) == 1:
            self._full_init(*args)
        elif len(args) == 10:
            self._short_init(*args)
        else:
            raise MagmaPlusInitException()

    def _full_init(self, line):
        parsed_line = parse_from_mgf(line)
        reference = scan = parsed_line[0]
        name = parsed_line[1]
        mass = float(parsed_line[2])
        mol_form = parsed_line[3]
        smiles = parsed_line[4]
        _ = parsed_line[5]
        inchi_key = parsed_line[6]
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            self._short_init(
                scan,
                mass,
                Chem.GetFormalCharge(mol),
                mol.GetNumHeavyAtoms(),
                Chem.MolToMolBlock(mol),
                inchi_key,
                mol_form,
                name,
                reference,
                Crippen.MolLogP(mol),
            )
        else:
            raise MagmaPlusInitException()

    def _short_init(
            self,
            scan,
            mass,
            charge,
            n_atoms,
            mol_block,
            inchi_key,
            mol_form,
            name,
            reference,
            logp,
    ):
        self.scan = scan
        self.mass = mass
        self.charge = charge
        self.n_atoms = n_atoms
        self.mol_block = mol_block
        self.inchi_key = inchi_key
        self.mol_form = mol_form
        self.name = name
        self.reference = reference
        self.logp = logp

    def to_tool_input(self, *args):
        pass
