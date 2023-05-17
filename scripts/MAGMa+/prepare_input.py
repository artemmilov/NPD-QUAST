import os
import sqlite3
import zlib

from rdkit import Chem
from rdkit.Chem import Crippen

with open(snakemake.output[1], 'w') as magma_ini:
    magma_ini.write(
        '''[magma job]
# Location of structure database to fetch candidate molecules to match against ms peak trees
structure_database.hmdb = {0}
chemical_engine = rdkit'''.format(snakemake.input[0])
    )

with open(snakemake.input[1]) as mgf:
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
    with open(snakemake.output[0], 'w') as mass_tree:
        mass_tree.write(record)
