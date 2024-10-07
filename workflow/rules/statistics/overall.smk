rule summary_table:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
        samples["path_to_sample"],
    output:
        "results/statistics/summary_table.html",
    log:
        "logs/summary_table.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/summary_table.py"


rule read_quality_histogram:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/read_quality_histogram.html",
    log:
        "logs/read_quality_histogram.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/read_quality_histogram.py"


rule passed_read_length_histogram:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/passed_read_length_histogram.html",
    log:
        "logs/passed_read_length_histogram.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/passed_read_length_histogram.py"


rule quality_per_base_position:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/quality_per_base_position.html",
    log:
        "logs/quality_per_base_position.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/quality_per_base_position.py"
