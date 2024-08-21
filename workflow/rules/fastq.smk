rule fastq:
    input:
        get_read_directory
    output:
        "results/fastq/{experiment_id}_{sample_id}.fastq.gz",
    log:
        "logs/fastq_{experiment_id}_{sample_id}.log",
    run:
        shell("""cp {input}/*.fastq.gz {output}""")
