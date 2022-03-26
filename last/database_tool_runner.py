import os.path
from subprocess import run

from rdkit import RDLogger

from last.NPDTools.npdtools_database import NpdToolsDatabase


def _prepare_input(challenge_file):
    new_filling = ''
    with open(challenge_file) as challenge:
        before_begin = True
        for line in challenge.readlines():
            if line == 'BEGIN IONS \n':
                new_filling += line
                before_begin = False
            elif line == 'END IONS \n':
                new_filling += line
            elif not before_begin:
                if '=' in line:
                    if line.split('=')[0].replace(' ', '') in ['PEPMASS', 'CHARGE', 'SCANS', 'MSLEVEL']:
                        new_filling += line
                elif '=' not in line:
                    new_filling += line
    with open(challenge_file, 'w') as challenge:
        challenge.write(new_filling)


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    candidates = range(70, 313)
    try:
        os.mkdir(
            os.path.join(
                '../sample',
                'temp',
            ),
        )
    except FileExistsError:
        pass
    database = NpdToolsDatabase(
        os.path.join(
            '../sample',
            'temp',
        ),
    )
    for candidate in candidates:
        path_to_database = os.path.join(
            '../sample',
            'CASMI2016_Cat2and3_Training_Candidates',
            f'Training-{candidate:03d}.csv',
        )
        path_to_challenge = os.path.join(
            '../sample',
            'CASMI2016_Cat2and3_Training_positive',
            f'Training-{candidate:03d}.mgf',
        )
        path_to_output = os.path.join(
            '../sample',
            'CASMI2016_Cat2and3_Training_Results',
            f'Training-{candidate:03d}',
        )
        _prepare_input(path_to_challenge)
        database.transform_from_casmi_data(path_to_database)
        run(
            [
                '/home/artem/Programming/bioinformatics/molDiscovery-2.6.0-beta-Linux/bin/dereplicator+.py',
                path_to_challenge,
                '--db-path',
                database.folder,
                '-o',
                path_to_output,
                '--pass-to-dereplicate',
                '--num_hits_to_report 10',
                '--min-score',
                '1',
                '--pm_thresh',
                '0.05'
            ],
        )
        database.clean()
    #os.rmdir(
    #    os.path.join(
    #        'files',
    #        'temp',
    #    ),
    #)


if __name__ == '__main__':
    main()
    # converter_from_NPDTools/output_sample/NPDTools_input      !!!
    # converter_from_NPDTools/merged_database      !!!
    # converter_from_NPDTools/test      !!!


#python MAGMa_plus.py read_ms_data

