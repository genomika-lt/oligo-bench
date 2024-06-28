rule process_fastq:
    input:
        directory=get_read_directory,
    output:
        "results/fastqs/{experiment_id}_{sample_id}.fastq.gz",
    log:
        "logs/fastqs_{experiment_id}_{sample_id}.log",
    run:
        input_dir = input.directory
        if not config["basecalling"]:
            shell(
                f"""
                cat {input_dir}/*.fastq.gz > {output}
            """
            )
        else:
            shell(
                f"""
                dorado basecaller sup@v5 \
                    --min-qscore 9 \
                    --no-trim \
                    --emit-fastq \
                    {input_dir} > {output}
            """
            )
