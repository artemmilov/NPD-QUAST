import os

specter = os.path.split(snakemake.input[0])[1]
challenge = os.path.split(os.path.split(os.path.split(snakemake.input[0])[0])[0])[1]

with open(
    os.path.join(
        snakemake.input[0],
        '0_{0}_FEATURE_1'.format(specter),
        'structure_candidates.tsv')
) as output:
    with open(snakemake.output[0], 'w') as tool_answers:
        for line in output.readlines()[1:]:
            answer_inchi_key = line.split('\t')[5]
            score = line.split('\t')[2]
            tool_answers.write(
                '{0}\t{1}\t{2}\t{3}\n'.format(
                    challenge,
                    specter,
                    answer_inchi_key,
                    str(-float(score)),
                ),
            )
