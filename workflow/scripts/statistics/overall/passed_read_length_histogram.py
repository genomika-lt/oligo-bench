"""Plots passed reads length histogram"""


import pysam
import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import parse_sam_records, snakemake_file_logger



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
        out = pysam.view('-o', 'out.sam', path_to_sample)
        records = [line.split() for line in out.split('\n')[:-1]]
        parsed_records = parse_sam_records(records)
        data = [len(record[10]) for record in parsed_records]
        figure.add_trace(go.Histogram(x=data,
                                      name=path_to_sample.split('/')[-1][7:-4]))

    figure.update_layout(title="Passed Read Length Histogram",
                         xaxis_title="Length",
                         yaxis_title="Number Of Reads",
                         legend_title="Samples")
    figure.update_layout(barmode='overlay')
    figure.update_traces(opacity=0.5)

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


passed_read_length_histogram(bam_files=snakemake.input, output_file=snakemake.output[0])
