rule dorado_basecaller:
    input:
        get_pod5_directory,
    output:
        "results/basecalled/{sample_id}.bam",
    log:
        "logs/basecall_{sample_id}.log",
    params:
        model=config["dorado_model"],
    resources:
        gpus=1,
    container:
        "docker://ontresearch/dorado:latest"
    shell:
        "dorado basecaller {params.model} --no-trim {input} > {output}"
