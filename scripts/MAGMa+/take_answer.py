import os
import sqlite3

specter = os.path.split(snakemake.output[0])[1].split('.')[0]
challenge = os.path.split(os.path.split(os.path.split(snakemake.output[0])[0])[0])[1]

with open(snakemake.output[0], 'w') as parsed_answers:
    conn = sqlite3.connect(snakemake.input[1])
    cur = conn.cursor()
    with open(snakemake.input[0]) as output:
        clean_output = []
        write = False
        for line in output.readlines():
            if (line == '') or (not line[0].isnumeric()):
                write = False
            if write:
                clean_output.append(line)
            if line == 'Candidate_Score Name Smiles\n':
                write = True

        for line in clean_output:
            answer_id = line.split(' ')[-2][1:-1]
            answer_inchi_key = list(
                cur.execute('SELECT * FROM molecules WHERE id = {0}'.format(answer_id))
            )[0][5]
            parsed_answers.write(
                '{0}\t{1}\t{2}\t{3}\n'.format(
                    challenge,
                    specter,
                    answer_inchi_key,
                    str(round(float(line.split(' ')[0]), 3)),
                ),
            )
