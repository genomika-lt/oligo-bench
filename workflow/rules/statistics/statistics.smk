include: "overall.smk"
include: "over_time.smk"


rule finalise_report:
    input:
        "results/statistics/summary_table.html",
        "results/statistics/mapping_table.html",
        "results/statistics/payload_table.html",
        "results/statistics/read_quality_histogram.html",
        "results/statistics/passed_read_length_histogram.html",
        "results/statistics/primers_number_histogram.html",
        # "results/statistics/quality_per_base_position.html",
        "results/statistics/pore_activity.html",
        "results/statistics/pore_scan.html",
    output:
        "results/report.html",
    log:
        "logs/finalise_report.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/finalise_report.py"
