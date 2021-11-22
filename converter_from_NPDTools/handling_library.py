import os
import shutil

from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit import RDLogger


from quast_mol import QuastMol, QuastMolInitException
from npdtools_database import NpdToolsDatabase, DatabaseInitException


def prepare_directories(
        database_folder,
        npdtools_input_file,
):
    try:
        shutil.rmtree(os.path.join(os.getcwd(), database_folder), ignore_errors=True)
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(os.getcwd(), database_folder))
    except OSError:
        print('Failed to make directory: {0}'.format(os.path.join(os.getcwd(), database_folder)))
        return False
    try:
        shutil.rmtree(os.path.join(os.getcwd(), npdtools_input_file), ignore_errors=True)
    except OSError:
        pass
    return True


def handle_library(
        library_file,
        database_folder,
        npdtools_input_file,
        answers_folder,
):
    database = NpdToolsDatabase(database_folder)
    with open(library_file, 'r', encoding='utf-8') as library, \
            open(npdtools_input_file, 'w') as npdtools_input, \
            open(os.path.join(answers_folder, 'true_answers.txt'), 'w') as true_answers:
        mols_data = library.read().split('\n\n')
        scan = 1
        for mol_data in mols_data:
            try:
                quast_mol = QuastMol(mol_data, scan)
            except QuastMolInitException as e:
                print(e)
                continue
            if database.add_mol(quast_mol):
                npdtools_input.write(
                    quast_mol.to_npdtools_input() + '\n\n',
                )
                true_answers.write(
                    '{0}\t{1}\n'.format(quast_mol.scan, quast_mol.inchi_key),
                )
                scan += 1


def print_dereplicator_answers(
        all_matches_file,
        dereplicator_answers_folder,
        database_folder,
):
    library_info_file = os.path.join(database_folder, 'library.info')
    with open(all_matches_file) as all_matches, \
            open(dereplicator_answers_folder, 'w') as dereplicator_answers, \
            open(library_info_file) as library_info:
        library_info_lines = library_info.readlines()
        for match in all_matches.readlines()[1:]:
            scan = match.split('\t')[2]
            answer_name = match.split('\t')[4]
            p_value = match.split('\t')[6]
            for library_info_line in library_info_lines:
                if library_info_line.split(' ')[1] == answer_name:
                    answer_address = os.path.join(database_folder, library_info_line.split(' ')[0])
                    answer_mol = Chem.MolFromMolFile(answer_address)
                    answer = inchi.MolToInchiKey(answer_mol).split('-')[0]
                    dereplicator_answers.write('{0}\t{1}\t{2}\n'.format(scan, answer, p_value))


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    while True:
        try:
            library_file, database_folder, npdtools_input_file, answers_folder = input().split()
            break
        except (ValueError, EOFError, SyntaxError):
            print('Incorrect input!')

    if not prepare_directories(database_folder, npdtools_input_file):
        return
    handle_library(
        library_file,
        database_folder,
        npdtools_input_file,
        answers_folder,
    )

    print('Handling library is done!')


if __name__ == '__main__':
    main()
    # input_sample/input_library.mgf output_sample/database output_sample/NPDTools_input/input.mgf test_main_input
    # dereplicator.py output_sample/NPDTools_input --db-path output_sample/database -o test
    # dereplicator.py output_sample/NPDTools_input/ --db-path pnpdatabase/ -o test
    # /home/artem/Programming/bioinformatics/molDiscovery-2.6.0-beta-Linux/bin/dereplicator.py
    # output_sample/NPDTools_input/ --db-path merged_database/ -o test --pass-to-dereplicate "--num_hits_to_report 10"
