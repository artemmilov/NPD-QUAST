import os.path
from subprocess import run

from rdkit import RDLogger

from converter_from_NPDTools.npdtools_database import NpdToolsDatabase


def main():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    candidates = range(59, 313)
    try:
        os.mkdir(
            os.path.join(
                'files',
                'temp',
            ),
        )
    except FileExistsError:
        pass
    database = NpdToolsDatabase(
        os.path.join(
            'files',
            'temp',
        ),
    )
    for candidate in candidates:
        path_to_database = os.path.join(
            'files',
            'CASMI2016_Cat2and3_Training_Candidates',
            f'Training-{candidate:03d}.csv',
        )
        path_to_challenge = os.path.join(
            'files',
            'CASMI2016_Cat2and3_Training_positive',
            f'Training-{candidate:03d}.mgf',
        )
        path_to_output = os.path.join(
            'files',
            'CASMI2016_Cat2and3_Training_Results',
            f'Training-{candidate:03d}',
        )
        database.transform_from_casmi_data(path_to_database)
        run(
            [
                '/home/artem/Programming/bioinformatics/molDiscovery-2.6.0-beta-Linux/bin/dereplicator.py',
                path_to_challenge,
                '--db-path',
                database.folder,
                '-o',
                path_to_output,
                '--pass-to-dereplicate',
                '--num_hits_to_report 10',
            ],
        )
        database.clean()
    os.rmdir(
        os.path.join(
            'files',
            'temp',
        ),
    )


if __name__ == '__main__':
    main()
    # converter_from_NPDTools/output_sample/NPDTools_input      !!!
    # converter_from_NPDTools/merged_database      !!!
    # converter_from_NPDTools/test      !!!
