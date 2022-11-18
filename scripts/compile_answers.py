with open(snakemake.output[0], 'w') as tool_answers:
    for ans_file in snakemake.input:
        print(snakemake.input)
        with open(ans_file) as cur:
            for line in cur.readlines():
                tool_answers.write(line)