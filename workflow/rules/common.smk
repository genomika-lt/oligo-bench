import os

import pandas as pd


samples = pd.read_csv("./config/experiments.csv").convert_dtypes()

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
    pod5_directory = os.path.join(path_to_sample.values[0], "./")
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


def get_fastq_files(wildcards):
    path_to_sample = samples.loc[
        samples["sample_id"] == wildcards.sample_id, "path_to_sample"
    ]
    path = path_to_sample.values[0]
    subdirectories = ["fastq", "fastq_pass", "fastq_fail"]
    fastq_files = []
    for subdir in subdirectories:
        fastq_directory = os.path.join(path, subdir)
        if os.path.exists(fastq_directory) and os.path.isdir(fastq_directory):
            for file in os.listdir(fastq_directory):
                if file.endswith(".fastq.gz") or file.endswith(".fastq"):
                    fastq_files.append(os.path.join(fastq_directory,file))

    for file in os.listdir(path):
        if file.endswith(".fastq.gz") or file.endswith(".fastq"):
            fastq_files.append(os.path.join(path,file))

    return fastq_files


def get_bam_files(wildcards):
    path_to_sample = samples.loc[
        samples["sample_id"] == wildcards.sample_id, "path_to_sample"
    ]
    path = path_to_sample.values[0]
    subdirectories = ["bam", "bam_pass", "bam_fail"]
    bam_files = []
    for subdir in subdirectories:
        bam_directory = os.path.join(path, subdir)
        if os.path.exists(bam_directory) and os.path.isdir(bam_directory):
            for file in os.listdir(bam_directory):
                if file.endswith(".bam"):
                    bam_files.append(os.path.join(bam_directory,file))

    for file in os.listdir(path):
        if file.endswith(".bam"):
            bam_files.append(os.path.join(path,file))

    return bam_files