import os

configfile: os.path.join('smk', 'config.yaml')


rule all:
    input:
       os.path.join(config['report_dir'], 'reports',
           config['report_name'], 'tool_answers.txt')


rule compile_answers:
    input:
        expand(os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge}.txt'), challenge=config['challenges'])
    output:
        os.path.join(config['report_dir'], 'reports',
            config['report_name'], 'tool_answers.txt')
    script:
        os.path.join('..', 'scripts', 'compile_answers.py')


rule take_answer_dereplicator_plus:
    input:
    output:


rule run_dereplicator_plus:
    input:
    output:


rule prepare_data_dereplicator_plus:
    input:
    output:
