""" Counts total number of reads inside SAM/BAM files """

import json
import os
from datetime import datetime

from pandas import DataFrame
from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_bam_file,
                           snakemake_file_logger)


# pylint: disable=too-many-locals
@snakemake_file_logger
def experiment_timestamps(path_to_samples, output_file):
    """
    Calculates start, end and duration of experiments and saves as csv file.
    :param list[str] path_to_samples: List with paths to experiments files.
    :param str output_file: Output file to save data.
    :rtype: None
    """


    data = {'Sample': [],
            'Start': [],
            'End': [],
            'Duration': []}

    for path in path_to_samples:
        path_to_json = [os.path.join(path, filename) for filename in os.listdir(path) if filename.endswith('.json')][0]

        with open(path_to_json, 'r', encoding='utf-8') as f:
            json_data = json.load(f)

        experiment_start_time = datetime.fromisoformat(json_data['protocol_run_info']['start_time'])

        experiment_end_time = datetime.fromisoformat(json_data['protocol_run_info']['end_time'])
        experiment_duration = (experiment_end_time - experiment_start_time).total_seconds()

        data['Sample'].append(path.split('/')[-2])
        data['Start'].append(experiment_start_time)
        data['End'].append(experiment_end_time)
        data['Duration'].append(experiment_duration)

    data_frame = DataFrame(data)
    data_frame.to_csv(output_file, index=False)


experiment_timestamps(path_to_samples=snakemake.input, output_file=snakemake.output[0])
