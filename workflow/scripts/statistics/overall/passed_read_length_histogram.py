"""Plots passed reads length histogram"""

import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger, parse_sam_bam_file


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
        for read in parse_sam_bam_file(path_to_sample):
            if len(read.query_sequence) not in data:
                data[len(read.query_sequence)] = 1
            else:
                data[len(read.query_sequence)] += 1
            number_of_reads += 1

        x = list(data.keys())
        x.sort()

        y = [data[key] for key in x]
        y_max = max(y)
        y_normalized = [i / y_max for i in y]

        filter_index_left = 0
        for i in range(len(x)):
            if y_normalized[i] > 5e-4:
                filter_index_left = i
                break

        filter_index_right = 0
        for i in range(len(x))[::-1]:
            if y_normalized[i] > 1e-2:
                filter_index_right = i
                break

        x_filtered = x[filter_index_left:filter_index_right]
        y_filtered = y_normalized[filter_index_left:filter_index_right]

        figure.add_trace(go.Scatter(y=y_filtered,
                                    x=x_filtered,
                                    name=path_to_sample.split('/')[-1][7:-4]))
    figure.update_xaxes(type="log")

    figure.update_layout(title="Passed Read Length Histogram",
                         xaxis_title="Length",
                         yaxis_title="Number Of Reads",
                         legend_title="Samples")

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


passed_read_length_histogram(bam_files=snakemake.input, output_file=snakemake.output[0])
