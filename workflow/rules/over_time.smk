rule gc_over_time:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/gc_over_time.html",
    log:
        "logs/gc_over_time.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/over_time/gc_over_time.py"


rule pore_activity:
    input:
        samples["path_to_sample"],
    output:
        "results/statistics/pore_activity.html",
    log:
        "logs/pore_activity.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/over_time/pore_activity.py"


rule pore_scan:
    input:
        samples["path_to_sample"],
    output:
        "results/statistics/pore_scan.html",
    log:
        "logs/pore_scan.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/over_time/pore_scan.py"


rule reads_over_time:
    input:
        samples["path_to_sample"],
    output:
        "results/statistics/reads_over_time.html",
    log:
        "logs/reads_over_time.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/over_time/reads_over_time.py"


rule quality_and_length_over_time:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/quality_and_length_over_time.html",
    log:
        "logs/quality_and_length_over_time.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/over_time/quality_and_length_over_time.py"


rule read_length_over_time_as_line:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/read_length_over_time_as_line.html",
    log:
        "logs/read_length_over_time_as_line.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/over_time/read_length_over_time_as_line.py"
