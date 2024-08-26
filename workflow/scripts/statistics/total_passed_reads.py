"""Counts total number of reads"""


import gzip
import logging

import plotly.graph_objects as go

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def count_total_passed_reads(fastq_files, output_file):
    """
    Counts total number of reads across samples
    :param str fastq_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    header = ['Parameter']
    values = ['Total reads']

    for file in fastq_files:
        with gzip.open(file, 'rt') as f:
            count = 0
            for _ in f.readlines():
                count += 1

        if count % 4:
            raise ValueError("Incorrect fastq file")
        count //= 4

        header.append(file[:-9])
        values.append(count)

    figure = go.Figure(data=[go.Table(header={'values': header},
                                      cells={'values': values})])

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html())


try:
    count_total_passed_reads(fastq_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
