rule test_multiqc_config:
    input:
        expand(
            "results/minknow/{experiment_id}_{sample_id}_throughput.csv",
            zip,
            experiment_id=samples["experiment_id"],
            sample_id=samples["sample_id"]
        ),
        config="config/multiqc_config.yaml",
    output:
        "results/qc/multiqc.html",
        "results/qc/multiqc.config_data.zip",
    params:
        extra="--verbose",
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.10.0/bio/multiqc"
