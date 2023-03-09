import os

configfile: os.path.join('smk', 'config.yaml')

rule all_dereplicator_plus:
    input:
       os.path.join(config['report_dir'], 'reports',
           config['report_name'], 'tool_answers.txt')

rule compile_answers_dereplicator_plus:
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
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge1}'),
        pre_result_folder_spectra=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'results', '{challenge1}', 'spectra')),
        pre_result_decoys=directory(os.path.join(
            config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge1}','decoys'))
    output:
        result=os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge1}.txt')
    script:
        os.path.join('..', 'scripts', 'Dereplicator+', 'take_answer.py')

rule run_dereplicator_plus:
    input:
        spectra=os.path.join(config['report_dir'], 'challenges', '{challenge2}', 'spectra'),
        decoys=os.path.join(config['report_dir'],'challenges','{challenge2}','decoys'),
        db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'databases' , '{challenge2}')
    output:
        res_spectra=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            config['run_tool'], 'results', '{challenge2}', 'spectra')),
        res_decoys=directory(os.path.join(
        config['report_dir'],'temp',config['report_name'],
            config['run_tool'],'results','{challenge2}','decoys'))
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
                '{output.res_spectra}'  # /path_to_result
            ] +
            config['options'] +
            [
                ';\n'
            ] +
            [
                os.path.join(
                    config['tool_dir'],
                    'bin',
                    'dereplicator+.py',
                ),
                '{input.decoys}',
                '--db-path',
                '{input.db}',
                '-o',
                '{output.res_decoys}'  # /path_to_result
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
