import os
import glob
import pandas as pd
from snakemake.utils import validate

samples = pd.read_csv(config["experiments"], dtype={"experiment_id": str, "sample_id": str, "run_dir": str})
validate(samples, schema="../schemas/experiments.schema.yaml")

wildcard_constraints:
    sample_id="|".join(samples["sample_id"]),
    experiment_id="|".join(samples["experiment_id"]),

def get_throughput_file(wildcards):
    try:
        run_dir = samples.loc[
            (samples["experiment_id"] == wildcards.experiment_id) & 
            (samples["sample_id"] == wildcards.sample_id), "run_dir"
        ].values[0]
    except IndexError:
        raise KeyError(f"No entry found for experiment_id={wildcards.experiment_id}, sample_id={wildcards.sample_id}") 
    throughput_pattern = os.path.join(run_dir, "throughput_*.csv")
    throughput_files = glob.glob(throughput_pattern)
    if throughput_files:
        return throughput_files[0]
    else:
        raise FileNotFoundError(f"No throughput file found for {wildcards.experiment_id}, {wildcards.sample_id}")

rule mod_throughput:
    input:
        get_throughput_file
    output:
        "results/minknow/{experiment_id}_{sample_id}_throughput.csv"
    log:
        "logs/mod_throughput_{experiment_id}_{sample_id}.log"
    params:
        samples=lambda wildcards: wildcards.experiment_id,
        id="Experiment Name"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/add_id_to_csv.py"
