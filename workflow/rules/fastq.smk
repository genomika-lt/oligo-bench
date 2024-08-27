rule fastq:
    input:
        get_fastq_directory,
    output:
        "results/basecalled/{experiment_id}_{sample_id}.fastq.gz",
    log:
        "logs/fastq_{experiment_id}_{sample_id}.log",
    run:
        shell("""cat {input}/*.fastq.gz > {output}""")
