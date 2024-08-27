rule basecall:
    input:
        get_pod5_directory,
    output:
        "results/basecalled/{experiment_id}_{sample_id}.bam",
    log:
        "logs/basecall_{experiment_id}_{sample_id}.log",
    run:
        shell(
            f"""dorado basecaller sup --min-qscore 9 \
                                        --no-trim \
                                        {input} > results/basecalled/{wildcards.experiment_id}_{wildcards.sample_id}.bam"""
        )
