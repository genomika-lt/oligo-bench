rule filter_basecalled:
    input:
        "results/basecalled/{sample_id}.bam",
    output:
        "results/basecalled/passed_{sample_id}.bam",
    log:
        "logs/filter_basecalled_{sample_id}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/filter_reads_quality.py"
