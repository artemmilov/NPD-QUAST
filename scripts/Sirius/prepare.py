import os

challenge_folder = os.path.split(os.path.split(snakemake.input[0])[0])[0]

with open(os.path.join(challenge_folder, 'database.csv')) as raw_database:
    with open(snakemake.output[0], 'w') as final_database:
        for line in raw_database.readlines():
            final_database.write(line)
