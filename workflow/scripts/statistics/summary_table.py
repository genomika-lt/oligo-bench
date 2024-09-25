"""Plots summary table and saves to html"""


import os
import json

from datetime import datetime

from workflow.scripts.utils import (parse_sam_records,
                                    file_logger,
                                    round_to_x_significant)

import pysam
import plotly.graph_objects as go

from snakemake.script import snakemake


def calculate_n50_for_numbers(list_of_lengths: list[int|float]):
    """
    Calculates N50 for list of lengths of sequences
    :param list[int] list_of_lengths: list with all lengths of sequences
    :return: Calculated N50
    :rtype: int
    """

    total_sum_half = sum(list_of_lengths) / 2
    for length in sorted(list_of_lengths, reverse=True):
        total_sum_half -= length
        if total_sum_half <= 0:
            return length

    return 0


@file_logger
def summary_table(bam_files, output_file):
    """
    Plots gc distribution over time in samples and saves to html
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    header_values = ['Sample ID',
                     'Total Reads, #', 'Passed Reads, %', 'Failed Reads, %',
                     'Passed Bases, #', 'N50', 'Passed Reads GC content, %', 'Experiment Time, h']
    body_values = [[], [], [], [], [], [], [], []]

    number_of_samples = len(bam_files) // 3
    files = [(bam_files[i + j * number_of_samples] for j in range(3)) for i in range(number_of_samples)]
    for file_all_reads, file_passed_reads, path_to_sample in files:
        out_all_reads = pysam.view('-o', 'out.sam', file_all_reads)
        total_reads_count = len(out_all_reads.split('\n')) - 1

        out_passed_reads = pysam.view('-o', 'out.sam', file_passed_reads)
        records_passed_reads = [record.split() for record in out_passed_reads.split('\n')[:-1]]
        parsed_passed_records = parse_sam_records(records_passed_reads)
        passed_reads_percentage = len(parsed_passed_records) / total_reads_count
        failed_reads_percentage = 1 - passed_reads_percentage
        passed_bases_count = sum(len(passed_record[9]) for passed_record in parsed_passed_records)
        passed_gc_percentage = sum(passed_record[9].count('G') + passed_record[9].count('C')
                                   for passed_record in parsed_passed_records) / passed_bases_count
        parsed_reads_lengths = [len(read[10]) for read in parsed_passed_records]
        n50 = calculate_n50_for_numbers(parsed_reads_lengths)

        json_file = [file for file in os.listdir(path_to_sample) if file.endswith('.json')][0]
        with open(os.path.join(path_to_sample, json_file), 'r') as f:
            json_data = json.load(f)

        experiment_start_time = datetime.strptime(json_data['protocol_run_info']['start_time'][:26],'%Y-%m-%dT%H:%M:%S.%f')
        experiment_end_time = datetime.strptime(json_data['protocol_run_info']['end_time'][:26],'%Y-%m-%dT%H:%M:%S.%f')
        experiment_duration_hours = (experiment_end_time - experiment_start_time).total_seconds() / 60 / 60

        body_values[0].append(file_all_reads.split('/')[-1].split('.')[0])
        body_values[1].append(total_reads_count)
        body_values[2].append(round_to_x_significant(passed_reads_percentage, 3))
        body_values[3].append(round_to_x_significant(failed_reads_percentage, 3))
        body_values[4].append(passed_bases_count)
        body_values[5].append(n50)
        body_values[6].append(round_to_x_significant(passed_gc_percentage, 3))
        body_values[7].append(round_to_x_significant(experiment_duration_hours, 3))

    figure = go.Figure(data=[go.Table(header={'values': header_values},
                                      cells={'values': body_values})])

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


summary_table(bam_files=snakemake.input, output_file=snakemake.output[0])
