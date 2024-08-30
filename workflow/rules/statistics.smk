rule count_total_passed_reads:
    input:
        expand(
            "results/basecalled/{experiment_id}_{sample_id}.bam",
            zip,
            experiment_id=samples["experiment_id"],
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
            "results/basecalled/{experiment_id}_{sample_id}.bam",
            zip,
            experiment_id=samples["experiment_id"],
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
            "results/basecalled/{experiment_id}_{sample_id}.bam",
            zip,
            experiment_id=samples["experiment_id"],
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
        samples['run_dir']
    output:
        "results/statistics/pore_activity.html",
    log:
        "logs/pore_activity.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/pore_activity.py"


rule finalise_report:
    input:
        "results/statistics/total_passed_reads.html",
        "results/statistics/calculate_n50.html",
        "results/statistics/gc_over_time.html",
        "results/statistics/pore_activity.html"
    output:
        "results/report.html",
    log:
        "logs/finalise_report.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/finalise_report.py"
