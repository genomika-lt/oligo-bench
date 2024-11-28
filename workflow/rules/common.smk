import os

import pandas as pd


samples = pd.read_csv(config["samples"]).convert_dtypes()

char = "/"  # workaround for snakemake linter false positive

samples.loc[:, "sample_id"] = [
    line.strip(char).split(char)[-2] for line in samples["path_to_sample"]
]


wildcard_constraints:
    sample_id="|".join(samples.loc[:, "sample_id"]),


def get_pod5_directory(wildcards):
    path_to_sample = samples.loc[
        samples["sample_id"] == wildcards.sample_id, "path_to_sample"
    ]
    pod5_directory = os.path.join(path_to_sample.values[0], "pod5")
    return pod5_directory


def get_reference(wildcards):
    path_to_reference = samples.loc[
        samples["sample_id"] == wildcards.sample_id, "path_to_reference"
    ].values[0]
    return path_to_reference


def get_forward_primer(wildcards):
    forward_primer = samples.loc[
        samples["sample_id"] == wildcards.sample_id, "forward_primer"
    ].values[0]
    return forward_primer


def get_reverse_primer(wildcards):
    reverse_primer = samples.loc[
        samples["sample_id"] == wildcards.sample_id, "reverse_primer"
    ].values[0]
    return reverse_primer
