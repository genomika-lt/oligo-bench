"""Plots reads quality histogram"""


import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import parse_sam_bam_file, snakemake_file_logger


@snakemake_file_logger
def read_quality_histogram(bam_files, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    figure = go.Figure()

    for path_to_sample in bam_files:
        data = {}
        number_of_reads = 0
        for read in parse_sam_bam_file(path_to_sample):
            quality_string = read.query_qualities
            mean_quality = sum(quality_string) // len(quality_string)

            if mean_quality not in data:
                data[mean_quality] = 1
            else:
                data[mean_quality] += 1

            number_of_reads += 1

        x = list(data.keys())
        x.sort()

        y = [data[key] for key in x]
        y_max = max(y)
        y_normalized = [i / y_max for i in y]

        figure.add_trace(go.Scatter(y=y_normalized,
                                    x=x,
                                    name=path_to_sample.split('/')[-1][:-4]))

    figure.update_layout(title="All Reads Mean Quality",
                         xaxis_title="Mean Quality",
                         yaxis_title="Number Of Reads",
                         legend_title="Samples")

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


read_quality_histogram(bam_files=snakemake.input, output_file=snakemake.output[0])
