"""Counts total number of reads"""


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


def count_total_passed_reads(bam_files, output_file):
    """
    Counts total number of reads across samples
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    header = ['Parameter']
    values = ['Total reads']

    for file in bam_files:
        out = pysam.view('-o', 'out.sam', file)
        count = len(out.split('\n')) - 1

        header.append(file.split('/')[-1])
        values.append(count)

    figure = go.Figure(data=[go.Table(header={'values': header},
                                      cells={'values': values})])

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


try:
    count_total_passed_reads(bam_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
