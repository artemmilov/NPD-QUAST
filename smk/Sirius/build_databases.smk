import os

configfile: os.path.join('smk', 'config.yaml')

rule all_sirius_build_databases:
    input:
        expand(os.path.join(config['report_dir'], 'temp', config['report_name'],
            'databases', '{challenge}.txt'), challenge=config['challenges'])


rule prepare_database_sirius_build_databases:
    input:
        raw_db=os.path.join(config['report_dir'], 'challenges', '{challenge}', 'database.csv')
    output:
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
                        'databases', '{challenge}.txt')
    script:
        os.path.join('..', '..', 'scripts', 'Sirius', 'prepare.py')
