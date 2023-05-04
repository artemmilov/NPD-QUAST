import os

configfile: os.path.join('smk','config.yaml')


rule all_magma_plus_run:
    input:
        os.path.join(config['report_dir'],'reports',
            config['report_name'],'tool_answers.txt')


def take_all_spectra(challenges):
    all_spectra = []
    for challenge in challenges:
        all_spectra += list(map(
            lambda s: os.path.join(challenge,'spectra',s.split('.')[0]),
            os.listdir(os.path.join(config['report_dir'],'challenges',
                challenge,'spectra'))
        ))
    return all_spectra


def form_data():
    challenges = []
    spectra = []
    for challenge in config['challenges']:
        for specter in os.listdir(os.path.join(config['report_dir'], 'challenges', challenge, 'spectra')):
            challenges.append(challenge)
            spectra.append(str(specter).split('.')[0])
    return [challenges, spectra]


rule compile_answers_magma_plus_run:
    input:
        expand(
            os.path.join(config['report_dir'],'temp',config['report_name'],
            'answers','{challenge}','spectra','{specter}.txt'),
            zip,
            challenge=form_data()[0],
            specter=form_data()[1]
        )
    output:
        os.path.join(config['report_dir'],'reports',
            config['report_name'],'tool_answers.txt')
    script:
        os.path.join('..','..','scripts','compile_answers.py')


rule take_answer_magma_plus_run:
    input:
        pre_result_folder=os.path.join(
            config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge}','spectra','{specter}','result.txt'),
        db=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{specter}.db')
    output:
        result=os.path.join(config['report_dir'],'temp',config['report_name'],
            'answers','{challenge}','spectra','{specter}.txt')
    script:
        os.path.join('..','..','scripts','MAGMa+','take_answer.py')


rule run_magma_plus_run:
    input:
        specter=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'inputs','{challenge}','spectra','{specter}.tree'),
        db=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'databases','{challenge}','spectra','{specter}.db'),
        job_ini=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge}','spectra','{specter}','magma_job.ini'),
        job_dir=os.path.join(config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge}','spectra','{specter}')
    output:
        res=os.path.join(
            config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge}','spectra','{specter}','result.txt')
    #    conda:
    #        '../envs/magma-plus-env.yml'
    shell:
        'cd {input.job_dir};\n' + \
        'export PATH={0}:$PATH;\n'.format(config['tool_dir']) + \
        #'mamba create -n kinoml --no-default-packages' + \
        #'mamba env update -n kinoml -f env.yml' + \
        #'eval "$(conda shell.bash hook)";\n' + \
        # 'conda env create -f ~/Programming/bioinformatics/NPD-QUAST/envs/magma-plus-env.yml || true;\n' + \
        'conda activate magma-plus-env-1;\n' + \
        'export MAGMAPLUS_CLASSIFIER_PATH={0};\n'.format(config['tool_dir']) + \
        'path_to_magma={0};\n'.format(config['tool_dir']) + \
        'python ' + config['tool_dir'] + '/MAGMa_plus.py' + \
        ' read_ms_data -i 1 -p 5 -q 0.001 -f mass_tree {input.specter} {input.db};\n' + \
        'python ' + config['tool_dir'] + '/MAGMa_plus.py' + \
        ' annotate -c 0 -d 0 -b 3 -w 1 -s hmdb {input.db};\n' + \
        'python ' + config['tool_dir'] + '/MAGMa_plus.py' + \
        ' export_result {input.db} > {output.res}'
