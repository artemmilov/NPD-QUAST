import os

configfile: os.path.join('smk', 'config.yaml')


def form_data():
    challenges = []
    spectra = []
    again_challenges = []
    decoys = []
    for challenge in config['challenges']:
        for specter in os.listdir(os.path.join(config['report_dir'], 'challenges', challenge, 'spectra')):
            challenges.append(challenge)
            spectra.append(str(specter).split('.')[0])
        for decoy in os.listdir(os.path.join(config['report_dir'],'challenges',challenge,'decoys')):
            again_challenges.append(challenge)
            decoys.append(str(decoy).split('.')[0])
    return [challenges, spectra, again_challenges, decoys]


rule all_sirius_run:
    input:
       os.path.join(config['report_dir'], 'reports',
           config['report_name'], 'tool_answers.txt')


# def take_all_spectra(challenges):
#     all_spectra = []
#     for challenge in challenges:
#         all_spectra += list(map(
#             lambda s: os.path.join(challenge, 'spectra', s.split('.')[0]),
#             os.listdir(os.path.join(config['report_dir'], 'challenges',
#                 challenge, 'spectra'))
#         ))
#     return all_spectra


rule compile_answers_sirius_run:
    input:
        [
            os.path.join(config['report_dir'],'temp',config['report_name'],
                'answers','{challenge}_{smt}_{specter}.txt').format(challenge=challenge,smt=smt, specter=specter)
            for challenge in config['challenges']
            for smt in os.listdir(os.path.join(config['report_dir'],'challenges',challenge))
            if os.path.isdir(os.path.join(config['report_dir'],'challenges',challenge, str(smt)))
            for specter in map(lambda s: s.split('.')[0],
            os.listdir(os.path.join(config['report_dir'],'challenges',challenge, str(smt))))
        ]
    output:
        os.path.join(config['report_dir'], 'reports',
            config['report_name'], 'tool_answers.txt')
    script:
         os.path.join('..', '..', 'scripts', 'compile_answers.py')


rule take_answer_sirius_run:
    input:
        pre_result_folder=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            'results', '{challenge}', '{smt}', '{specter}'))
    output:
        result=os.path.join(config['report_dir'], 'temp', config['report_name'],
            'answers', '{challenge}_{smt}_{specter}.txt')
    script:
        os.path.join('..', 'scripts', 'Sirius', 'take_answer.py')


def get_challenge_by_specter(specter):
    _temp = os.path.split(os.path.split(specter)[0])[0]
    if _temp[-1] == '/':
        return _temp[:-1]
    return temp


rule run_sirius_run:
    input:
        specter=os.path.join(config['report_dir'], 'challenges', '{challenge}', '{smt}', '{specter}' + '.mgf'),
        dirty_db=os.path.join(config['report_dir'], 'temp', config['report_name'],
            'databases' , '{specter}') + '.txt'
    output:
        res=directory(os.path.join(
            config['report_dir'], 'temp', config['report_name'],
            'results', '{challenge}', '{smt}', '{specter}'))
    shell:
        'export PATH=\"{0}:$PATH\";'.format(os.path.join(config['tool_dir'], 'bin')) + config['a']['d']['e'] +\
        'sirius -i ' + '{input.dirty_db}' + ' {0} custom-db {1};'.format(
            config['options'][0], config['options'][1]) + \
        'sirius -i {input.specter} -o {output.res} ' + '{0} formula {1} structure {2}'.format(
            config['options'][2], config['options'][3], config['options'][4])

# export PATH="/home/artem/Programming/bioinformatics/sirius-5.5.7-linux64-headless/sirius/bin:$PATH"
# sirius -i /home/artem/Programming/bioinformatics/NPD-QUAST-test/snakemake_test/new2_unexecuted/temp/database.txt  custom-db --name cur_database
# sirius -i /home/artem/Programming/bioinformatics/NPD-QUAST-test/snakemake_test/new2_unexecuted/temp/spectra/Challenge-082.mgf -o /home/artem/Programming/bioinformatics/NPD-QUAST-test/snakemake_test/new2_unexecuted/temp/tool/cur_results/Challenge-082  formula -c 10 structure --database cur_database

