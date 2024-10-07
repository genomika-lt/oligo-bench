"""Plots mean base quality for each position"""


import pysam
import plotly.express as px
from pandas import DataFrame

from workflow.scripts.utils import parse_sam_records, snakemake_file_logger, from_char_to_phred

from snakemake.script import snakemake


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
        out = pysam.view('-o', 'out.sam', path_to_sample)
        records = [line.split() for line in out.split('\n')[:-1]]
        parsed_records = parse_sam_records(records)

        qualities = [from_char_to_phred(record[10]) for record in parsed_records]
        sample_data = [0 for i in range(100)]

        for record in qualities:
            for percentage in range(100):
                start = int(len(record) / 100 * percentage)
                end = int(len(record) / 100 * (percentage + 1))
                if start == end:
                    continue
                sample_data[percentage] += sum(record[start:end]) / (end - start)

        data_for_plotting[path_to_sample.split('/')[-1][7:-4]] = [i / len(qualities) for i in sample_data]


    df = DataFrame(data_for_plotting)
    figure = px.line(df, x=df.index, y=df.columns)

    figure.update_layout(xaxis_title='Percentage Position',
                         yaxis_title='Quality',
                         title='Mean Quality per percentage position',
                         legend_title='Samples')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


quality_per_base_position(bam_files=snakemake.input, output_file=snakemake.output[0])
