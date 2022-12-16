import os
import sqlite3
import zlib

from rdkit import Chem
from rdkit.Chem import Crippen

with open(snakemake.output[2], 'w') as magma_ini:
    magma_ini.write(
        '''[magma job]
# Location of structure database to fetch candidate molecules to match against ms peak trees
structure_database.hmdb = {0}
chemical_engine = rdkit'''.format(snakemake.output[0])
    )

with open(snakemake.input[0]) as mgf:
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
    with open(snakemake.output[1], 'w') as mass_tree:
        mass_tree.write(record)


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
        raise RuntimeError('Wrong smiles: {0}'.format(smiles))
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


conn = sqlite3.connect(snakemake.output[0])
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
except Exception:
    pass  # ????????????????????????????

challenge_folder = os.path.split(os.path.split(snakemake.input[0])[0])[0]
with open(os.path.join(challenge_folder, 'database.csv')) as raw_database:
    for line in raw_database.readlines()[1:]:
        if line != '':
            try:
                c.execute(
                    '''INSERT INTO molecules (id, mim, charge, natoms, molblock, inchikey,
                     molform, name, reference, logp) VALUES (?,?,?,?,?,?,?,?,?,?)''',
                    _parse_line_for_magma(line),
                )
            except Exception:
                pass
conn.commit()

c.execute('PRAGMA temp_store = 2')
c.execute('''CREATE INDEX idx_cover ON molecules (charge, mim, natoms, reference, 
molform, inchikey, name, molblock, logp)''')
conn.commit()
