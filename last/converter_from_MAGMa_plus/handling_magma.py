import glob
import os
from subprocess import run

from database_creator import create_database
from mass_tree_creator import create_mass_tree


def parse_input():
    while True:
        try:
            work_folders = input().split()
            return work_folders
        except (ValueError, EOFError, SyntaxError):
            print('Incorrect input!')


def prepare_folder(folder):
    work_folder = os.path.join(folder, 'work')
    files = glob.glob(work_folder)
    for f in files:
        os.remove(f)
    with open(os.path.join(work_folder, 'magma_job.ini'), 'w') as ini:
        ini.write('[magma job]\n\n')
        ini.write('# Location of structure database to fetch ')
        ini.write('candidate molecules to match against ms peak trees\n')
        ini.write(
            'structure_database.hmdb = {0}\n\n'.format(
                'database.db',
            ),
        )
        ini.write('chemical_engine = rdkit')


def write_true_answers(folder):
    pass


def transform_data(folder):
    pass


def run_magma(folder):
    pass


def parse_output(folder):
    pass


def write_tool_answers(output):
    pass


def handle_reporting(folder):
    pass


def main():
    folders = parse_input()
    for folder in folders:
        prepare_folder(folder)
        write_true_answers(folder)
        transform_data(folder)
        run_magma(folder)
        output = parse_output(folder)
        write_tool_answers(output)
        handle_reporting(folder)


def main():
    try:
        os.mkdir('handling_magma_temp')
    except FileExistsError:
        os.rmdir('handling_magma_temp')
        os.mkdir('handling_magma_temp')
    database_file, spectra_file = handle_input()
    sql_database_file = os.path.join(
        'handling_magma_temp',
        database_file.split('.')[0] + '.db',
    )
    mass_tree_file = os.path.join(
        'handling_magma_temp',
        spectra_file.split('.')[0] + '.tree',
    )
    create_database(database_file, sql_database_file)
    create_mass_tree(spectra_file, mass_tree_file)
    path_to_magma = os.path.relpath(
        os.curdir,
        '/home/artem/Programming/bioinformatics/MAGMa-plus',
    )
    run(
        [
            'export',
            'MAGMAPLUS_CLASSIFIER_PATH={0}'.format(path_to_magma),
        ],
    )
    run(
        [
            'python',
            os.path.join(
                path_to_magma,
                'MAGMa_plus.py',
            ),
            'read_ms_data',
            '-i',
            '1',
            '-p',
            '5',
            '-q',
            '0.001',
            '-f',
            'mass_tree',
            mass_tree_file,
            'output.db',
        ],
    )
    run(
        os.path.join(
            'handling_magma_temp',
            database_file.split('.')[0],
        ),
    )


if __name__ == '__main__':
    main()
