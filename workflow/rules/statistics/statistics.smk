include: "overall.smk"
include: "over_time.smk"


rule finalise_report:
    input:
        f"results/statistics/summary_table.html",
        f"results/statistics/mapping_table.html",
        f"results/statistics/payload_table.html",
        f"results/statistics/read_quality_histogram.html",
        f"results/statistics/passed_read_length_histogram.html",
        f"results/statistics/primers_number_histogram.html",
        # "results/statistics/quality_per_base_position.html",
        f"results/statistics/pore_activity.html",
        f"results/statistics/pore_scan.html",
    output:
        f"{config["settings"]["save_path"]}/report.html",
    log:
        "logs/finalise_report.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/finalise_report.py"
