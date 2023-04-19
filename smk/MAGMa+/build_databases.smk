import os

configfile: os.path.join('smk', 'config.yaml')

rule all_magma_plus_build_databases:
    input:
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge}.db'), challenge=config['challenges']),
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{specter}.tree')),
        expand(os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{specter}','magma_job.ini'))


rule prepare_data_magma_plus_build_databases:
    input:
        raw_spectrer = os.path.join(config['report_dir'],'challenges','{specter}.mgf')
    output:
        db=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{specter}.db'),
        spectrer=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{specter}.tree'),
        job_ini=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{specter}','magma_job.ini')
    script:
        os.path.join('..', '..', 'scripts', 'MAGMa+', 'prepare.py')
