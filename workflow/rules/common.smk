import os
import glob
import pandas as pd
from snakemake.utils import validate

samples = pd.read_csv(
    config["experiments"],
    dtype={"experiment_id": str, "sample_id": str, "run_dir": str},
)

validate(samples, schema="../schemas/experiments.schema.yaml")


wildcard_constraints:
    sample_id="|".join(samples["sample_id"]),
    experiment_id="|".join(samples["experiment_id"]),


def get_throughput_file(wildcards):
    try:
        run_dir = samples.loc[
            (samples["experiment_id"] == wildcards.experiment_id)
            & (samples["sample_id"] == wildcards.sample_id),
            "run_dir",
        ].values[0]
    except IndexError:
        raise KeyError(
            f"No entry found for experiment_id={wildcards.experiment_id}, sample_id={wildcards.sample_id}"
        )
    throughput_pattern = os.path.join(run_dir, "throughput_*.csv")
    throughput_files = glob.glob(throughput_pattern)
    if throughput_files:
        return throughput_files[0]
    else:
        raise FileNotFoundError(
            f"No throughput file found for {wildcards.experiment_id}, {wildcards.sample_id}"
        )


def get_pore_activity_file(wildcards):
    try:
        run_dir = samples.loc[
            (samples["experiment_id"] == wildcards.experiment_id)
            & (samples["sample_id"] == wildcards.sample_id),
            "run_dir",
        ].values[0]
    except IndexError:
        raise KeyError(
            f"No entry found for experiment_id={wildcards.experiment_id}, sample_id={wildcards.sample_id}"
        )
    pore_pattern = os.path.join(run_dir, "pore_activity_*.csv")
    pore_files = glob.glob(pore_pattern)
    if pore_files:
        return pore_files[0]
    else:
        raise FileNotFoundError(
            f"No pore activity file found for {wildcards.experiment_id}, {wildcards.sample_id}"
        )


def get_read_directory(wildcards):
    try:
        run_dir = samples.loc[
            (samples["experiment_id"] == wildcards.experiment_id)
            & (samples["sample_id"] == wildcards.sample_id),
            "run_dir",
        ].values[0]
    except IndexError:
        raise KeyError(
            f"No entry found for experiment_id={wildcards.experiment_id}, sample_id={wildcards.sample_id}"
        )
    if config["basecalling"]:
        input_dir = os.path.join(run_dir, "pod5_pass")
    else:
        input_dir = os.path.join(run_dir, "fastq_pass")
    if not os.path.exists(input_dir):
        raise ValueError(
            f"Directory {input_dir} does not exist for {wildcards.experiment_id}/{wildcards.sample_id}"
        )
    return input_dir
