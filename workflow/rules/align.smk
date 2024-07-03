rule minimap2:
    input:
        target=config["reference"],
        query="results/fastqs/{experiment_id}_{sample_id}.fastq.gz",
    output:
        "results/aligned/{experiment_id}_{sample_id}_aln.sorted.bam",
    log:
        "logs/minimap2_{experiment_id}_{sample_id}.log",
    params:
        extra="-ax sr -L --MD -t8 -Y -k10 -w5",
        sorting="queryname",
        sort_extra="",
    threads: 3
    wrapper:
        "v3.13.3/bio/minimap2/aligner"

rule samtools_stats:
    input:
        bam="results/aligned/{experiment_id}_{sample_id}_aln.sorted.bam",
    output:
        "results/aligned/{experiment_id}_{sample_id}.stats",
    params:
        extra="",
    log:
        "logs/samtools_stats_{experiment_id}_{sample_id}.log",
    wrapper:
        "v3.13.3/bio/samtools/stats"
