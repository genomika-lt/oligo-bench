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


rule finalise_report:
    input:
        "results/statistics/total_passed_reads.html",
        "results/statistics/calculate_n50.html"
    output:
        "results/report.html",
    log:
        "logs/finalise_report.log",
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/statistics/finalise_report.py"
