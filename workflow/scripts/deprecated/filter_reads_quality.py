"""Plots summary table and saves to html"""


import json
import os
import logging

from datetime import datetime

import pysam
import plotly.graph_objects as go

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


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
            meta_key, meta_type, meta_value = meta_data.split(':')
            if meta_type == 'i':
                meta[meta_key] = int(meta_value)
            elif meta_type == 'f':
                meta[meta_key] = float(meta_value)
            elif meta_type == 'Z':
                meta[meta_key] = meta_value
            else:
                raise ValueError(f"Unexpected Data Type: {meta_type}")

        parsed_record = record[:11] + [meta]
        parsed_records.append(parsed_record)

    return parsed_records


def filter_basecalled(bam_file, output_file):
    """
    Plots gc distribution over time in samples and saves to html
    :param str bam_file: Folder with fastq reads
    :param str output_file: Output file to save data to
    :rtype: None
    """
    print(bam_file)
    data = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
    filtered = pysam.AlignmentFile(output_file, "wb", template=data, check_sq=False)
    for read in data:
        for tag in read.tags:
            if tag[0] == 'qs':
                if tag[1] > 9:
                    filtered.write(read)

try:
    filter_basecalled(bam_file=snakemake.input[0], output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
