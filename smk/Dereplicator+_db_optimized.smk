import os

configfile: os.path.join('smk', 'config.yaml')

rule all_dereplicator_plus_db_optimized:
    input:
       os.path.join(config['report_dir'], 'reports',
           config['report_name'], 'tool_answers.txt')

rule compile_answers_dereplicator_plus_db_optimized:
    input:
        expand(os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge}'), challenge=config['challenges'])
    output:
        os.path.join(config['report_dir'], 'reports',
            config['report_name'], 'tool_answers.txt')
    script:
        os.path.join('..', 'scripts', 'compile_answers.py')

rule get_challenge_answers_db_optimized:
    input:
        expand(os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge}', ''), challenge=config['challenges'])
    output:
        os.path.join(config['report_dir'], 'temp', config['report_name'], 'answers', '{challenge}')
    script:
        os.path.join('..', 'scripts', 'compile_answers.py')









rule take_answer_dereplicator_plus_db_optimized:
    input:
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge1}'),
        pre_result_folder=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'results', '{challenge1}'))
    output:
        result=os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge1}.txt')
    script:
        os.path.join('..', 'scripts', 'Dereplicator+', 'take_answer.py')

rule run_dereplicator_plus_db_optimized:
    input:
        spectra=os.path.join(config['report_dir'], 'challenges', '{challenge2}', 'spectra'),
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge2}')
    output:
        res=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'results', '{challenge2}'))
    shell:
        ' '.join(
            [
                os.path.join(
                    config['tool_dir'],
                    'bin',
                    'dereplicator+.py',
                ),
                '{input.spectra}',
                '--db-path',
                '{input.db}',
                '-o',
                '{output.res}'  # /path_to_result
            ] +
            config['options']
        )


rule prepare_data_dereplicator_plus_db_optimized:
    input:
        raw_db=os.path.join(config['report_dir'], 'challenges', '{challenge3}', 'database.csv')
    output:
        db=directory(os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge3}'))
    script:
        os.path.join('..', 'scripts', 'Dereplicator+', 'prepare.py')
