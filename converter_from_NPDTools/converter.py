import shutil

from rdkit import Chem
from rdkit.Chem import Descriptors
import os


def main():
    try:
        gnps_library_address, npdtools_database_address, npdtools_input_address = input().split()
    except (ValueError, EOFError):
        print('Incorrect input!')
        return

    try:
        shutil.rmtree(os.getcwd() + '/' + npdtools_database_address, ignore_errors=True)
    except OSError:
        pass
    try:
        os.mkdir(os.getcwd() + '/' + npdtools_database_address)
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + npdtools_database_address)
        return
    try:
        os.mkdir(os.getcwd() + '/' + npdtools_database_address + '/mols')
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + npdtools_database_address + '/mols')
        return
    try:
        shutil.rmtree(os.getcwd() + '/' + npdtools_input_address, ignore_errors=True)
    except OSError:
        pass
    try:
        os.mkdir(os.getcwd() + '/' + npdtools_input_address)
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + npdtools_input_address)
        return

    with open(gnps_library_address + '.mgf', 'r', encoding='utf-8') as file_gnps_library, \
            open(npdtools_database_address + '/library.info', 'w') as file_npdtools_database_library, \
            open(npdtools_database_address + '/smiles.info', 'w') as file_smiles, \
            open(npdtools_input_address + '.mgf', 'w') as file_npdtools_input_address:
        mols_data = file_gnps_library.read().split('\n\n')
        for i, mol_data in enumerate(mols_data):
            mol_data_for_input = []
            smiles, filename, name, mass, num_acids = None, None, None, None, None
            filename = f'compound_{i:05d}'
            num_acids = 1000

            for j, s in enumerate(filter(lambda _s: _s != '', mol_data.splitlines())):
                if '=' in s:
                    s_split = s.split('=', maxsplit=1)
                    left, right = s_split[0], s_split[1]
                    if (left == 'SMILES') and ('.' not in right):
                        smiles = right
                    elif (left == 'NAME') and (left != ''):
                        name = right.replace(' ', '_')
                    elif left == 'PEPMASS':
                        mass = float(right)
                    if j <= 3:
                        mol_data_for_input.append(s)
                else:
                    mol_data_for_input.append(s)
            if (filename is None) or (name is None) or (mass is None) or (num_acids is None):
                continue

            try:
                m = Chem.MolFromSmiles(smiles)
                if m is None:
                    continue
                m = Chem.AddHs(m)
                Chem.MolToMolFile(m, npdtools_database_address + '/mols/' + filename + '.mol', forceV3000=True)
                file_npdtools_database_library.write(
                    'mols/{}.mol {} {} {} DB\n'.format(filename, name, Descriptors.ExactMolWt(m), num_acids))
                file_smiles.write(smiles + '\n')

                for s in mol_data_for_input:
                    file_npdtools_input_address.write(s + '\n')
                file_npdtools_input_address.write('\n\n')
            except Exception as e:
                print(e)

    print('Success!')


if __name__ == '__main__':
    main()
