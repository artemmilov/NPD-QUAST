import os

configfile: os.path.join('smk', 'config.yaml')

def form_data():
    challenges = []
    spectra = []
    for challenge in config['challenges']:
        for specter in os.listdir(os.path.join(config['report_dir'], 'challenges', challenge, 'spectra')):
            challenges.append(challenge)
            spectra.append(str(specter).split('.')[0])
    return [challenges, spectra]

rule all_magma_plus_build_databases:
    input:
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge}.db'), challenge=config['challenges']),
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{challenge1}','spectra','{specter}.tree'),
            challenge1=form_data()[0],
            specter=form_data()[1]
        ),
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge1}','spectra','{specter}','magma_job.ini'),
            challenge1=form_data()[0],
            specter=form_data()[1])


rule prepare_data_magma_plus_build_databases:
    input:
        raw_specter = os.path.join(config['report_dir'],'challenges','{challenge1}','spectra','{specter}.mgf')
    output:
        db=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge1}.db'),
        specter=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{challenge1}','spectra','{specter}.tree'),
        job_ini=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge1}','spectra','{specter}','magma_job.ini')
    script:
        os.path.join('..', '..', 'scripts', 'MAGMa+', 'prepare.py')
