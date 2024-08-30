rule throughput:
    input:
        get_throughput_file,
    output:
        "results/minknow/{experiment_id}_{sample_id}_throughput.csv",
    log:
        "logs/throughput_{experiment_id}_{sample_id}.log",
    params:
        samples=lambda wildcards: wildcards.experiment_id,
        id="Experiment Name",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/add_id_to_csv.py"


rule pore_activity:
    input:
        get_pore_activity_file,
    output:
        "results/minknow/{experiment_id}_{sample_id}_pore_act.csv",
        "results/minknow/{experiment_id}_{sample_id}_pore_seq.csv",
        "results/minknow/{experiment_id}_{sample_id}_pore_ratio.csv",
    log:
        "logs/pore_{experiment_id}_{sample_id}.log",
    params:
        samples=lambda wildcards: wildcards.experiment_id,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/parse_pore_activity.py"