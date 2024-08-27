rule minimap2:
    input:
        target=config["reference"],
        query="results/basecalled/{experiment_id}_{sample_id}.bam",
    output:
        "results/aligned/{experiment_id}_{sample_id}_aln.sorted.bam",
    log:
        "logs/minimap2_{experiment_id}_{sample_id}.log",
    params:
        extra="-ax sr -L --MD -t8 -Y -k10 -w5 -m10",
        sorting="coordinate",
        sort_extra="",
    threads: 4
    wrapper:
        "v3.13.3/bio/minimap2/aligner"


rule samtools_stats:
    input:
        bam="results/aligned/{experiment_id}_{sample_id}_aln.sorted.bam",
    output:
        "results/aligned/stats/{experiment_id}_{sample_id}.stats",
    params:
        extra="",
    log:
        "logs/samtools_stats_{experiment_id}_{sample_id}.log",
    wrapper:
        "v3.13.3/bio/samtools/stats"


rule parse_samtools_stats:
    input:
        "results/aligned/stats/{experiment_id}_{sample_id}.stats",
    output:
        "results/aligned/stats/{experiment_id}_{sample_id}_MAPQ.csv",
        "results/aligned/stats/{experiment_id}_{sample_id}_INS.csv",
        "results/aligned/stats/{experiment_id}_{sample_id}_DEL.csv",
        "results/aligned/stats/{experiment_id}_{sample_id}_COV.csv",
    log:
        "logs/parse_stats_{experiment_id}_{sample_id}.log",
    params:
        samples=lambda wildcards: wildcards.experiment_id,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/parse_samtools_stats.py"


rule count_seqs_bam:
    input:
        bam="results/aligned/{experiment_id}_{sample_id}_aln.sorted.bam",
        ref=config["reference"],
    output:
        out="results/aligned/stats/{experiment_id}_{sample_id}_ref_dist.csv",
    params:
        samples=lambda wildcards: wildcards.experiment_id,
    log:
        "logs/count_sequences_{experiment_id}_{sample_id}.log",
    threads: 4
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/count_seqs_bam.py"


rule count_mism_bam:
    input:
        bam="results/aligned/{experiment_id}_{sample_id}_aln.sorted.bam",
    output:
        out="results/aligned/stats/{experiment_id}_{sample_id}_nm_counts.csv",
    log:
        "logs/count_mism_{experiment_id}_{sample_id}.log",
    threads: 4
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/count_mism_bam.py"
