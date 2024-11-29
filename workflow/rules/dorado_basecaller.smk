rule dorado_basecaller:
    input:
        get_pod5_directory,
    output:
        "results/basecalled/{sample_id}.bam",
    log:
        "logs/basecall_{sample_id}.log",
    params:
        model=config["basecalling"]["dorado_model"]["value"],
    resources:
        gpus=1,
    conda:
        "../envs/pysam.yaml"
    shell:
        "dorado basecaller {params.model} -r {input} >> {output}"
