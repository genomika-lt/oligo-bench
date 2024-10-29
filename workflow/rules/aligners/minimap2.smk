rule minimap2:
    input:
        target=get_reference,
        query="results/basecalled/passed_{sample_id}.bam",
    output:
        "results/aligned/{sample_id}.bam",
    log:
        "logs/minimap2_{sample_id}.log",
    params:
        extra="-ax sr -L --MD -Y -k10 -w5 -m10",
        sorting="coordinate",
        sort_extra="",
    threads: 4
    wrapper:
        "v3.13.3/bio/minimap2/aligner"


rule samtools_stats:
    input:
        bam="results/aligned/{sample_id}.bam",
    output:
        "results/aligned/stats/{sample_id}.stats",
    params:
        extra="",
    log:
        "logs/samtools_stats_{sample_id}.log",
    wrapper:
        "v3.13.3/bio/samtools/stats"


rule parse_samtools_stats:
    input:
        "results/aligned/stats/{sample_id}.stats",
    output:
        "results/aligned/stats/{sample_id}_MAPQ.csv",
        "results/aligned/stats/{sample_id}_DEL.csv",
        "results/aligned/stats/{sample_id}_INS.csv",
        "results/aligned/stats/{sample_id}_COV.csv",
    log:
        "logs/parse_stats_{sample_id}.log",
    params:
        samples=lambda wildcards: wildcards.sample_id,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/parse_samtools_stats.py"


rule count_seqs_bam:
    input:
        bam="results/aligned/{sample_id}.bam",
        ref=get_reference,
    output:
        out="results/aligned/stats/{sample_id}_ref_dist.csv",
    params:
        samples=lambda wildcards: wildcards.sample_id,
    log:
        "logs/count_sequences_{sample_id}.log",
    threads: 4
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/count_seqs_bam.py"


rule count_mism_bam:
    input:
        bam="results/aligned/{sample_id}.bam",
    output:
        out="results/aligned/stats/{sample_id}_nm_counts.csv",
    log:
        "logs/count_mism_{sample_id}.log",
    threads: 4
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/count_mism_bam.py"
