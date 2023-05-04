import os
from rdkit import Chem

id_to_inchi = {}
with open(os.path.join(snakemake.input[0], 'library.smiles')) as smiles:
    for i, line in enumerate(smiles.readlines()):
        if line != '':
            m = Chem.MolFromSmiles(line)
            if m is not None:
                id_to_inchi[i] = Chem.MolToInchiKey(m).split('-')[0]
            else:
                id_to_inchi[i] = 'ERROR'

with open(snakemake.output[0], 'w') as tool_answers:
    with open(os.path.join(snakemake.input[1], 'all_matches.tsv')) as output:
        for line in output.readlines()[1:]:
            # QUESTION!!!QUESTION!!!QUESTION!!!QUESTION!!!QUESTION!!!
            challenge_name = os.path.split(os.path.split(os.path.split(line.split('\t')[0])[0])[0])[1]
            # QUESTION!!!QUESTION!!!QUESTION!!!QUESTION!!!QUESTION!!!
            spectra_or_decoy = os.path.split(os.path.split(line.split('\t')[0])[0])[1]
            answer_scan = int(line.split('\t')[2])
            answer_id = int(line.split('\t')[3])
            answer_inchi_key = id_to_inchi[answer_id]
            specter = os.path.split(line.split('\t')[0])[-1].split('.')[0]
            score = line.split('\t')[5]
            tool_answers.write(
                '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    challenge_name,
                    spectra_or_decoy,
                    specter,
                    answer_scan,
                    answer_inchi_key,
                    str(-float(score)),
                ),
            )
