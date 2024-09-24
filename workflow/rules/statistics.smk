rule count_total_passed_reads:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            zip,
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total_passed_reads.html",
    log:
        "logs/count_total_passed_reads.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/total_passed_reads.py"


rule calculate_n50:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/calculate_n50.html",
    log:
        "logs/calculate_n50.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/N50.py"


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
        "../scripts/statistics/gc_over_time.py"


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
        "../scripts/statistics/pore_activity.py"


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
        "../scripts/statistics/pore_scan.py"


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
        "../scripts/statistics/reads_over_time.py"


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
        "../scripts/statistics/quality_and_length_over_time.py"


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
        expand(
            "results/aligned/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
        samples["path_to_sample"],
    output:
        "results/statistics/summary_table.html",
    log:
        "logs/summary_table.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/summary_table.py"


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
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/read_quality_histogram.py"


rule finalise_report:
    input:
        "results/statistics/summary_table.html",
        "results/statistics/read_quality_histogram.html",
        "results/statistics/quality_and_length_over_time.html",
        "results/statistics/total_passed_reads.html",
        "results/statistics/pore_activity.html",
        "results/statistics/calculate_n50.html",
        "results/statistics/gc_over_time.html",
        "results/statistics/pore_scan.html",
        "results/statistics/reads_over_time.html",
    output:
        "results/report.html",
    log:
        "logs/finalise_report.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/finalise_report.py"
