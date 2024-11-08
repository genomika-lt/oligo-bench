"""Plots passed reads length histogram"""

import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import parse_sam_records, snakemake_file_logger, parse_bam_file



@snakemake_file_logger
def passed_read_length_histogram(bam_files, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    # pylint: disable=duplicate-code
    figure = go.Figure()

    for path_to_sample in bam_files:
        data = {}
        number_of_reads = 0
        for read in parse_bam_file(path_to_sample):
            if len(read.qual) not in data:
                data[len(read.qual)] = 1
            else:
                data[len(read.qual)] += 1
            number_of_reads += 1

        x = list(data.keys())
        x.sort()

        y = [data[key] for key in x]
        y_max = max(y)
        y_normalized = [i / y_max for i in y]

        figure.add_trace(go.Scatter(y=y_normalized,
                                    x=x,
                                    name=path_to_sample.split('/')[-1][7:-4]))
    figure.update_xaxes(type="log")

    figure.update_layout(title="Passed Read Length Histogram",
                         xaxis_title="Length",
                         yaxis_title="Number Of Reads",
                         legend_title="Samples")
    figure.update_layout(barmode='overlay')
    figure.update_traces(opacity=0.65)

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


passed_read_length_histogram(bam_files=snakemake.input, output_file=snakemake.output[0])
