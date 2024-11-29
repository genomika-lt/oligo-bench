"""Plots mean base quality for each position"""


import pysam
import plotly.express as px

from pandas import DataFrame
from pysam.bcftools import query

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger, parse_sam_bam_file



# pylint: disable=too-many-locals
@snakemake_file_logger
def quality_per_base_position(bam_files, output_file):
    """
    Plots mean base quality for each position in percents
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    data_for_plotting = {}

    for path_to_sample in bam_files:
        sample_data = [0 for _ in range(100)]

        read_counter = 0
        for read in parse_sam_bam_file(path_to_sample):
            read_counter += 1
            percent_read_length = len(read.query_qualities)

            for percentage in range(100):
                start = round(percent_read_length * percentage)
                end = round(percent_read_length * (percentage + 1))
                if start == end:
                    continue
                sample_data[percentage] += sum(read.query_qualities[start:end]) / (end - start)

        data_for_plotting[path_to_sample.split('/')[-1][7:-4]] = [i / read_counter
                                                                  for i in sample_data]

    df = DataFrame(data_for_plotting)
    figure = px.line(df, x=df.index, y=df.columns)

    figure.update_layout(xaxis_title='Percentage Position',
                         yaxis_title='Quality',
                         title='Mean Quality per percentage position',
                         legend_title='Samples')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


quality_per_base_position(bam_files=snakemake.input, output_file=snakemake.output[0])
