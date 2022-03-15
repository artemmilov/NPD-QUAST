import os
import shutil

from rdkit import Chem

from converter_from_Sirius.sirius_tool import SiriusTool
from general import parse_from_mgf
from npd_quast_folder import NPDQuastFolder


def main():
    path_to_spectres = \
        '../../data/NPD-QUAST/CASMI2016_Cat2and3_Training_positive'
    path_to_databases = \
        '../../data/NPD-QUAST/CASMI2016_Cat2and3_Challenge_Candidates'
    path_to_quast_folder = \
        '../sample/large'

    for spectra in os.listdir(path_to_spectres):
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
                try:
                    smiles = Chem.MolToSmiles(
                        Chem.AddHs(Chem.MolFromSmiles(parse_from_mgf(dataline)[4])),
                    )
                except Exception:
                    continue
                database.write(smiles + '\n')


def add_sirius_report():
    path_to_npd_quast_folder = '../sample/large'
    npd_quast_folder = NPDQuastFolder(path_to_npd_quast_folder)
    sirius_tool = SiriusTool()
    npd_quast_folder.make_tool_report(sirius_tool)


if __name__ == '__main__':
    main()
    #add_sirius_report()
