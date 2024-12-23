rule provide_bam:
    params:
        fastq_files = get_fastq_files,
        bam_files= get_bam_files
    output:
        "results/basecalled/{sample_id}.bam"
    log:
        "logs/convert_fastq_bam_{sample_id}.log"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/fastq_bam_convert.py"

