import os

from converter_from_MAGMa_plus.magma_tool import MagmaTool
from mass_tree_creator import create_mass_tree
from database_creator import create_database
from npd_quast_folder import NPDQuastFolder


def main():
    path_to_spectres = \
        '../files/CASMI2016_Cat2and3_Training_positive'
    path_to_databases = \
        '../files/CASMI2016_Cat2and3_Challenge_Candidates'
    path_to_quast_folder = \
        '../files/casmi_quast'

    for spectra in os.listdir(path_to_spectres)[1:2]:
        name = spectra.split('.')[0]
        os.mkdir(
            os.path.join(path_to_quast_folder, 'challenges', name),
        )
        os.mkdir(
            os.path.join(path_to_quast_folder, 'challenges', name, 'spectres'),
        )
        create_mass_tree(
            os.path.join(path_to_spectres, name + '.mgf'),
            os.path.join(path_to_quast_folder, 'challenges', name, 'spectres', name + '.tree'),
        )
        create_database(
            os.path.join(path_to_databases, name + '.csv'),
            os.path.join(path_to_quast_folder, 'challenges', name, name + '.db'),
        )


def add_magma_plus_report():
    path_to_npd_quast_folder = '../files/casmi_quast'
    npd_quast_folder = NPDQuastFolder(path_to_npd_quast_folder)
    tool = MagmaTool()
    npd_quast_folder.make_tool_report(tool)


if __name__ == '__main__':
    main()
    add_magma_plus_report()
