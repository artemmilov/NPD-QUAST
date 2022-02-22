import os
import shutil

from converter_from_Sirius.sirius_tool import SiriusTool
from npd_quast_folder import NPDQuastFolder
from general import parse_from_mgf


def main():
    path_to_spectres = \
        '../files/CASMI2016_Cat2and3_Training_positive'
    path_to_databases = \
        '../files/CASMI2016_Cat2and3_Challenge_Candidates'
    path_to_quast_folder = \
        '../files/casmi_quast'

    for spectra in os.listdir(path_to_spectres)[:1]:
        name = spectra.split('.')[0]
        if os.path.isdir(os.path.join(path_to_quast_folder, 'challenges', name)):
            shutil.rmtree(os.path.join(path_to_quast_folder, 'challenges', name))
        os.mkdir(os.path.join(path_to_quast_folder, 'challenges', name))
        os.mkdir(os.path.join(path_to_quast_folder, 'challenges', name, 'spectres'))
        shutil.copyfile(
            os.path.join(path_to_spectres, name + '.mgf'),
            os.path.join(path_to_quast_folder, 'challenges', name, 'spectres', name + '.mgf'),
        )
        with open(os.path.join(path_to_databases, name + '.csv')) as raw_database, \
                open(
                    os.path.join(
                        path_to_quast_folder,
                        'challenges',
                        name,
                        name + '.txt',
                    ),
                    'w',
                ) as database:
            for dataline in raw_database.readlines()[1:]:
                if dataline == '':
                    continue
                smiles = parse_from_mgf(dataline)[4]
                database.write(smiles + '\n')


def add_sirius_report():
    path_to_npd_quast_folder = '../files/casmi_quast'
    npd_quast_folder = NPDQuastFolder(path_to_npd_quast_folder)
    sirius_tool = SiriusTool()
    npd_quast_folder.make_tool_report(sirius_tool)


if __name__ == '__main__':
    main()
    add_sirius_report()
