import os

configfile: os.path.join('smk', 'config.yaml')

rule all_dereplicator_plus:
    input:
       os.path.join(config['report_dir'], 'reports',
           config['report_name'], 'tool_answers.txt')


rule compile_answers_dereplicator_plus:
    input:
        [
            os.path.join(config['report_dir'], 'temp', config['report_name'],
                'answers', '{challenge}_{smt}.txt').format(challenge=challenge, smt=smt)
            for challenge in config['challenges']
            for smt in os.listdir(os.path.join(config['report_dir'], 'challenges', challenge))
            if os.path.isdir(os.path.join(config['report_dir'], 'challenges', challenge, str(smt)))
         ]
    output:
        os.path.join(config['report_dir'], 'reports',
            config['report_name'], 'tool_answers.txt')
    script:
        os.path.join('..', 'scripts', 'compile_answers.py')


rule take_answer_dereplicator_plus:
    input:
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge}'),
        pre_result_folder=os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'results', '{challenge}', '{smt}')
    output:
        result=os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge}_{smt}.txt')
    script:
        os.path.join('..', 'scripts', 'Dereplicator+', 'take_answer.py')


rule run_dereplicator_plus:
    input:
        inp=os.path.join(config['report_dir'], 'challenges', '{challenge2}', '{smt}'),
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge2}')
    output:
        res_spectra=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'results', '{challenge2}', '{smt}'))
    shell:
        ' '.join(
            [
                os.path.join(
                    config['tool_dir'],
                    'bin',
                    'dereplicator+.py',
                ),
                '{input.inp}',
                '--db-path',
                '{input.db}',
                '-o',
                '{output.res_spectra}'  # /path_to_result
            ] +
            config['options']
        )


rule prepare_data_dereplicator_plus:
    input:
        raw_db=os.path.join(config['report_dir'], 'challenges', '{challenge3}', 'database.csv')
    output:
        db=directory(os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge3}'))
    script:
        os.path.join('..', 'scripts', 'Dereplicator+', 'prepare.py')
