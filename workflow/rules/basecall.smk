rule dorado_basecaller:
    input:
        get_pod5_directory,
    output:
        "results/basecalled/{sample_id}.bam",
    log:
        "logs/basecall_{sample_id}.log",
    params:
        model=config['dorado_model']
    shell:
        "dorado basecaller {params.model} --no-trim {input} > {output}"
