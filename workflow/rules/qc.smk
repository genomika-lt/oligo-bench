rule fastqc:
    input:
        "results/fastqs/{experiment_id}_{sample_id}.fastq.gz",
    output:
        html="results/qc/fastqc/{experiment_id}_{sample_id}.html",
        zip="results/qc/fastqc/{experiment_id}_{sample_id}_fastqc.zip",
    log:
        "logs/fastqc_{experiment_id}_{sample_id}.log",
    threads: 1
    resources:
        mem_mb=1024,
    wrapper:
        "v3.10.0/bio/fastqc"


rule multiqc:
    input:
        expand(
            "results/minknow/{experiment_id}_{sample_id}_throughput.csv",
            zip,
            experiment_id=samples["experiment_id"],
            sample_id=samples["sample_id"],
        ),
        expand(
            "results/qc/fastqc/{experiment_id}_{sample_id}_fastqc.zip",
            zip,
            experiment_id=samples["experiment_id"],
            sample_id=samples["sample_id"],
        ),
        expand(
            "results/aligned/stats/{experiment_id}_{sample_id}.stats",
            zip,
            experiment_id=samples["experiment_id"],
            sample_id=samples["sample_id"],
        ),
        config="config/multiqc_config.yaml",
    output:
        "results/qc/multiqc.html",
    params:
        extra="--verbose",
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.10.0/bio/multiqc"
