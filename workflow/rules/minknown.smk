rule throughput:
    input:
        get_throughput_file,
    output:
        "results/minknow/{experiment_id}_{sample_id}_throughput.csv",
    log:
        "logs/mod_throughput_{experiment_id}_{sample_id}.log",
    params:
        samples=lambda wildcards: wildcards.experiment_id,
        id="Experiment Name",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/add_id_to_csv.py"
