rule basecall:
    input:
        get_pod5_directory
    output:
        "results/fastq/{experiment_id}_{sample_id}.fastq.gz",
    log:
        "logs/basecall_{experiment_id}_{sample_id}.log",
    run:
        shell(f"""dorado basecaller sup --min-qscore 9 \
                                        --no-trim \
                                        --emit-fastq \
                                        {input} > results/fastq/{wildcards.experiment_id}_{wildcards.sample_id}.fastq &&
                  gzip results/fastq/{wildcards.experiment_id}_{wildcards.sample_id}.fastq""")
