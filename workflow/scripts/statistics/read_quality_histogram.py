"""Plots reads quality histogram"""

import logging
import sys

import pysam
import plotly.graph_objects as go

from workflow.scripts.utils.parse import parse_sam_records

from snakemake.script import snakemake

# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def read_quality_histogram(bam_files, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    figure = go.Figure()

    for path_to_sample in bam_files:
        out = pysam.view('-o', 'out.sam', path_to_sample)
        records = [line.split() for line in out.split('\n')[:-1]]
        parsed_records = parse_sam_records(records)
        data = [record[11]['qs'] for record in parsed_records]
        figure.add_trace(go.Histogram(x=data,
                                      name=path_to_sample.split('/')[-1][:-4]))

    figure.update_layout(title="Read Quality Histogram",
                         xaxis_title="Quality",
                         yaxis_title="Number Of Reads",
                         legend_title="Samples")
    figure.update_layout(barmode='overlay')
    figure.update_traces(opacity=0.5)

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


try:
    read_quality_histogram(bam_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
