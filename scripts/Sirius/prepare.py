with open(snakemake.input[0]) as raw_database:
    with open(snakemake.output[0], 'w') as final_database:
        for line in raw_database.readlines():
            final_database.write(line)
