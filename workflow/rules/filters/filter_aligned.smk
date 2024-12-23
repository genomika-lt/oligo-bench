rule filter_aligned_reads:
    input:
        "results/aligned/{sample_id}.bam",
    output:
        "results/aligned/passed_{sample_id}.bam",
    log:
        "logs/{sample_id}/filter_aligned_reads.log",
    conda:
        "../../envs/pysam.yaml"
    shell:
        "samtools view -h {input} --min-MQ 10 --excl-flags 260 >> {output}"
