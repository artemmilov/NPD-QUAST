import os

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi


def print_dereplicator_answers(npdtools_output_folder,
                               database_folder,
                               answers_folder):
    with open(os.path.join(npdtools_output_folder, 'all_matches.tsv')) as all_matches, \
            open(os.path.join(database_folder, 'library.info')) as library_info, \
            open(os.path.join(answers_folder, 'tool_answers.txt'), 'w') as dereplicator_answers:
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
    while True:
        try:
            npdtools_output_folder, database_folder, answers_folder = input().split()
            break
        except (ValueError, EOFError, SyntaxError):
            print('Incorrect input!')

    print_dereplicator_answers(npdtools_output_folder,
                               database_folder,
                               answers_folder)

    print('Handling NPDTools is done!')


if __name__ == '__main__':
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    main()  # test output_sample/database test_main_input
