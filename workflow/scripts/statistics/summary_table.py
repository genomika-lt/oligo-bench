"""Plots summary table and saves to html"""


import json
import os
import logging


from math import log10, floor
from datetime import datetime

import pysam
import pandas as pd
import plotly.graph_objects as go

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def round_to_x_significant(number, x):
    return round(number, x - 1 - floor(log10(abs(number))))

def parse_bam_records(records):
    """
    Parse bam records.
    :param list[list[str]] records: Records from BAM file.
    :return: Parsed records.
    """

    parsed_records = []
    for record in records:
        meta = {}
        for meta_data in record[11:]:
            key_index = meta_data.index(':')
            type_index = meta_data.index(':', key_index + 1)
            meta_key, meta_type, meta_value = (meta_data[:key_index],
                                               meta_data[key_index + 1:type_index],
                                               meta_data[type_index + 1:])
            if meta_type == 'i':
                meta[meta_key] = int(meta_value)
            elif meta_type == 'f':
                meta[meta_key] = float(meta_value)
            elif meta_type == 'Z':
                meta[meta_key] = meta_value
            else:
                logger.warning("Unexpected Data Type: %s", meta_type)
                meta[meta_key] = meta_value

        parsed_record = record[:11] + [meta]
        parsed_records.append(parsed_record)

    return parsed_records


def summary_table(bam_files, output_file):
    """
    Plots gc distribution over time in samples and saves to html
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    header_values = ['Sample ID',
                     'Total Reads, #', 'Passed Reads, %', 'Failed Reads, %', 'Mapped Reads, %',
                     'Passed Bases, #', 'Mapped Bases, #', 'Mapped Bases, %',
                     'Passed Reads GC content, %', 'Experiment Time, h']
    body_values = [[], [], [], [], [], [], [], [], [], []]

    number_of_samples = len(bam_files) // 4
    files = [(bam_files[i + j * number_of_samples] for j in range(4)) for i in range(number_of_samples)]

    for file_all_reads, file_passed_reads, file_mapped_reads, path_to_sample in files:
        out_all_reads = pysam.view('-o', 'out.sam', file_all_reads)
        total_reads_count = len(out_all_reads.split('\n')) - 1

        out_passed_reads = pysam.view('-o', 'out.sam', file_passed_reads)
        records_passed_reads = [record.split() for record in out_passed_reads.split('\n')[:-1]]
        parsed_passed_records = parse_bam_records(records_passed_reads)
        passed_reads_percentage = len(parsed_passed_records) / total_reads_count
        failed_reads_percentage = 1 - passed_reads_percentage
        passed_bases_count = sum(len(passed_record[9]) for passed_record in parsed_passed_records)
        passed_gc_percentage = sum(passed_record[9].count('G') + passed_record[9].count('C')
                                   for passed_record in parsed_passed_records) / passed_bases_count

        out_mapped_reads = pysam.view('-o', 'out.sam', file_mapped_reads)
        records_mapped_reads = [record.split() for record in out_mapped_reads.split('\n')[:-1]]
        parsed_mapped_records = parse_bam_records(records_mapped_reads)
        mapped_reads_percentage = len(parsed_mapped_records) / len(parsed_passed_records)
        mapped_bases_count = sum(len(mapped_record[9]) for mapped_record in parsed_mapped_records)
        mapped_bases_percentage = mapped_bases_count / passed_bases_count

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
        body_values[4].append(round_to_x_significant(mapped_reads_percentage, 3))
        body_values[5].append(passed_bases_count)
        body_values[6].append(mapped_bases_count)
        body_values[7].append(round_to_x_significant(mapped_bases_percentage, 3))
        body_values[8].append(round_to_x_significant(passed_gc_percentage, 3))
        body_values[9].append(round_to_x_significant(experiment_duration_hours, 3))

    figure = go.Figure(data=[go.Table(header={'values': header_values},
                                      cells={'values': body_values})])

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


try:
    summary_table(bam_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
