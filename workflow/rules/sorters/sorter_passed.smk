rule sort_passed_reads:
    input:
        "results/basecalled/passed_{sample_id}.bam",
    output:
        "results/sorted/sorted_{sample_id}.bam",
    log:
        "logs/{sample_id}/sort_passed_reads.log",
    conda:
        "../../envs/pysam.yaml"
    shell:
        "samtools sort {input} -t st -T temp_sorted_{input} -o {output}"
