import shutil

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import inchi
import os


def prepare_directories(npdtools_database_address,
                        npdtools_input_address,
                        true_answers_address):
    try:
        shutil.rmtree(os.path.join(os.getcwd(), npdtools_database_address), ignore_errors=True)
    except OSError:
        pass
    try:
        os.mkdir(os.getcwd() + '/' + npdtools_database_address)
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + npdtools_database_address)
        return False
    try:
        os.mkdir(os.getcwd() + '/' + npdtools_database_address + '/mols')
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + npdtools_database_address + '/mols')
        return False
    try:
        shutil.rmtree(os.getcwd() + '/' + npdtools_input_address, ignore_errors=True)
    except OSError:
        pass
    try:
        os.mkdir(os.getcwd() + '/' + npdtools_input_address)
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + npdtools_input_address)
        return False
    try:
        shutil.rmtree(os.getcwd() + '/' + true_answers_address, ignore_errors=True)
    except OSError:
        pass
    try:
        os.mkdir(os.getcwd() + '/' + true_answers_address)
    except OSError:
        print('Failed to make directory: ' + os.getcwd() + '/' + true_answers_address)
        return False
    return True


def handle_mols_data(gnps_library_address,
                     npdtools_database_address,
                     npdtools_input_address,
                     true_answers_address):
    with open(gnps_library_address + '.mgf', 'r', encoding='utf-8') as file_gnps_library, \
            open(npdtools_database_address + '/library.info', 'w') as file_npdtools_database_library, \
            open(npdtools_database_address + '/smiles.info', 'w') as file_smiles, \
            open(npdtools_input_address + '.mgf', 'w') as file_npdtools_input, \
            open(true_answers_address + '/true_answers.txt', 'w') as file_true_answers:
        mols_data = file_gnps_library.read().split('\n\n')
        scan = 1
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
                    file_npdtools_input.write(s + '\n')
                    if s == 'BEGIN IONS':
                        file_npdtools_input.write('SCANS={0}\n'.format(scan))
                file_npdtools_input.write('\n\n')

                file_true_answers.write('{}\t{}\n'.format(scan, inchi.MolToInchiKey(m).split('-')[0]))

                scan += 1
            except Exception as e:
                print(e)


def print_dereplicator_answers(all_matches_address,
                               dereplicator_answers_address,
                               database_address):
    library_info_address = os.path.join(database_address, 'library.info')
    with open(all_matches_address) as all_matches, \
            open(dereplicator_answers_address, 'w') as dereplicator_answers, \
            open(library_info_address) as library_info:
        library_info_lines = library_info.readlines()
        for match in all_matches.readlines()[1:]:
            scan = match.split('\t')[2]
            answer_name = match.split('\t')[4]
            p_value = match.split('\t')[6]
            for library_info_line in library_info_lines:
                if library_info_line.split(' ')[1] == answer_name:
                    answer_address = os.path.join(database_address, library_info_line.split(' ')[0])
                    answer_mol = Chem.MolFromMolFile(answer_address)
                    answer = inchi.MolToInchiKey(answer_mol).split('-')[0]
                    dereplicator_answers.write('{0}\t{1}\t{2}\n'.format(scan, answer, p_value))


def main():
    try:
        gnps_library_address, npdtools_database_address, npdtools_input_address, true_answers_address = input().split()
    except (ValueError, EOFError):
        print('Incorrect input!')
        return

    if not prepare_directories(npdtools_database_address, npdtools_input_address, true_answers_address):
        return

    handle_mols_data(gnps_library_address,
                     npdtools_database_address,
                     npdtools_input_address,
                     true_answers_address)

    print('Success!')

    print_dereplicator_answers('test/all_matches.tsv',
                               'tool_answers.txt',
                               'test/library_info.tsv')


if __name__ == '__main__':
    #print_dereplicator_answers('test/all_matches.tsv',
    #                           'tool_answers.txt',
    #                           'output_sample/database')
    main()  # input_sample/input_library output_sample/database output_sample/NPDTools_input output_sample/true_answers
    # dereplicator.py output_sample/NPDTools_input/ --db-path output_sample/database/ -o test
