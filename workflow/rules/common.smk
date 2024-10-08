import os
import glob

import pandas as pd

from snakemake.utils import validate


samples = pd.read_csv(config["samples"]).convert_dtypes()
samples.loc[:, "sample_id"] = [
    line.split("/")[-2] for line in samples["path_to_sample"]
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
