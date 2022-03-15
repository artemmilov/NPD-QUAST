import os
import shutil

from converter_from_NPDTools.dereplicator_plus_tool import DereplicatorPlusTool
from converter_from_NPDTools.dereplicator_tool import DereplicatorTool
from npd_quast_folder import NPDQuastFolder


def main():
    path_to_spectres = \
        '../sample/CASMI2016_Cat2and3_Training_positive'
    path_to_databases = \
        '../sample/CASMI2016_Cat2and3_Challenge_Candidates'
    path_to_quast_folder = \
        '../sample/large'

    for spectra in os.listdir(path_to_spectres)[:1]:
        name = spectra.split('.')[0]
        if os.path.isdir(os.path.join(path_to_quast_folder, 'challenges', name)):
            shutil.rmtree(os.path.join(path_to_quast_folder, 'challenges', name))
        os.mkdir(os.path.join(path_to_quast_folder, 'challenges', name))
        os.mkdir(os.path.join(path_to_quast_folder, 'challenges', name, 'spectres'))
        with open(
                os.path.join(path_to_spectres, spectra)
        ) as raw_spectra, open(
            os.path.join(path_to_quast_folder, 'challenges', name, 'spectres', spectra),
            'w',
        ) as converted_spectra:
            converted_spectra.write('BEGIN IONS\n')
            converted_spectra.write('MSLEVEL=2\n')
            for line in raw_spectra.readlines():
                if line == '':
                    continue
                if str.isnumeric(line[0]):
                    converted_spectra.write(line)
                else:
                    if '=' in line:
                        if line.split('=')[0] in ['SCANS', 'PEPMASS', 'CHARGE']:
                            converted_spectra.write(line)
            converted_spectra.write('END IONS\n')
        shutil.copyfile(
            os.path.join(path_to_databases, name + '.csv'),
            os.path.join(path_to_quast_folder, 'challenges', name, name + '.csv'),
        )


def add_dereplicator_plus_report():
    path_to_npd_quast_folder = '../sample/large'
    npd_quast_folder = NPDQuastFolder(path_to_npd_quast_folder)
    dereplicator_plus_tool = DereplicatorPlusTool()
    dereplicator_tool = DereplicatorTool()
    npd_quast_folder.make_tool_report(dereplicator_plus_tool, default_ranks=[1, 2, 3])
    npd_quast_folder.make_tool_report(dereplicator_tool, default_ranks=[4, 5, 6])


if __name__ == '__main__':
    # main()
    add_dereplicator_plus_report()
