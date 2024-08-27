"""Calculates N50 for each sample and saves to html table"""


import logging

import pysam
import plotly.graph_objects as go

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def calculate_n50_for_numbers(list_of_lengths):
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

def calculate_n50_for_bam_files(bam_files, output_file):
    """
    Calculates N50 for each sample and saves to html table
    :param str bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    header = ['Parameter']
    values = ['N50']

    for file in bam_files:
        out = pysam.view('-o', 'out.sam', file)
        list_of_lengths = [len(line.split('\t')[9]) for line in out.split('\n')[:-1]]

        n50 = calculate_n50_for_numbers(list_of_lengths)

        header.append(file.split('/')[-1])
        values.append(n50)

    figure = go.Figure(data=[go.Table(header={'values': header},
                                      cells={'values': values})])

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False))


try:
    calculate_n50_for_bam_files(bam_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
