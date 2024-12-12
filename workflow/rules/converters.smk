rule passed_bam_to_fastq:
    input:
        "results/basecalled/passed_{sample_id}.bam",
    output:
        "results/basecalled/passed_{sample_id}.fastq",
    log:
        "logs/passed_bam_to_fastq_{sample_id}.log",
    conda:
        "../envs/pysam.yaml"
    shell:
        "samtools fastq {input} > {output}"
