rule fastq:
    input:
        get_pod5_directory
    output:
        "results/fastq/{experiment_id}_{sample_id}.fastq.gz",
    log:
        "logs/fastq_{experiment_id}_{sample_id}.log",
    run:
        shell(f"""dorado basecaller sup@v5 --min-qscore 9 \
                                           --no-trim \
                                           --emit-fastq \
                                           {input} > {output}""")
