import os

configfile: os.path.join('smk', 'config.yaml')

def form_data():
    challenges = []
    smts = []
    spectra = []
    for challenge in config['challenges']:
        for specter in os.listdir(os.path.join(config['report_dir'], 'challenges', challenge, 'spectra')):
            challenges.append(challenge)
            smts.append('spectra')
            spectra.append(str(specter).split('.')[0])
        if os.path.exists(os.path.join(config['report_dir'], 'challenges', challenge, 'decoys')):
            for specter in os.listdir(os.path.join(config['report_dir'], 'challenges', challenge, 'decoys')):
                challenges.append(challenge)
                smts.append('decoys')
                spectra.append(str(specter).split('.')[0])
    return [challenges, smts, spectra]

rule all_magma_plus_build_databases:
    input:
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge}.db'), challenge=config['challenges']),
        expand(
            os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{challenge}','{smt}','{specter}.tree'),
            zip,
            challenge=form_data()[0],
            smt=form_data()[1],
            specter=form_data()[2]
        ),
        expand(
            os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge}','{smt}','{specter}','magma_job.ini'),
            zip,
            challenge=form_data()[0],
            smt=form_data()[1],
            specter=form_data()[2]
        )


rule prepare_input_magma_plus_build_databases:
    input:
        db = os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge}.db'),
        raw_specter = os.path.join(config['report_dir'],'challenges','{challenge}','{smt}','{specter}.mgf')
    output:
        specter=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{challenge}','{smt}','{specter}.tree'),
        job_ini=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge}','{smt}','{specter}','magma_job.ini')
    script:
        os.path.join('..', '..', 'scripts', 'MAGMa+', 'prepare_input.py')


rule prepare_database_magma_plus_build_databases:
    input:
        raw_db = os.path.join(config['report_dir'],'challenges','{challenge}','database.csv')
    output:
        db=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge}.db')
    script:
        os.path.join('..', '..', 'scripts', 'MAGMa+', 'prepare_database.py')
