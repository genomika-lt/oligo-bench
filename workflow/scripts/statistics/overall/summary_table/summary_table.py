"""Plots summary table for samples from precalculated csv files and saves to html. """


import os

import pysam
import pandas as pd
import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_records,
                           snakemake_file_logger,
                           round_to_x_significant,
                           integer_to_human_readable)

def ration_to_percentage_string(ration: float, digits_after_point=2) -> str:
    return f'{100 * round_to_x_significant(ration, 2 + digits_after_point):.{digits_after_point}f}%'

# pylint: disable=too-many-locals
@snakemake_file_logger
def summary_table(csv_files, output_file):
    """
    Plots summary table for samples from precalculated csv files and saves to html.
    :param dict[str, str] csv_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    total_reads_number = pd.read_csv(csv_files['total_reads'])
    passed_reads_number = pd.read_csv(csv_files['passed_reads'])
    total_bases_number = pd.read_csv(csv_files['total_bases'])
    passed_bases_number = pd.read_csv(csv_files['passed_bases'])
    gc_passed_bases_number = pd.read_csv(csv_files['gc_passed_bases'])
    timestamps = pd.read_csv(csv_files['timestamps'])

    bases_alphabet = ('b', 'Kb', 'Mb', 'Gb', 'Tb', 'Pb')
    passed_reads_ratio = (passed_reads_number.loc[:, 'Passed Reads']
                          / total_reads_number.loc[:, 'Total Reads'])
    gc_content_ratio = (gc_passed_bases_number.loc[:, 'GC Bases']
                         / passed_bases_number.loc[:, 'Passed Bases'])

    header_values = ['Sample ID',
                     'Total Reads',
                     'Total Bases',
                     'Passed Reads',
                     'Passed Bases',
                     'Passed GC',
                     'Experiment Duration']

    body_values = [total_reads_number.loc[:, 'Sample'],
                   total_reads_number.loc[:, 'Total Reads'].apply(integer_to_human_readable),
                   total_bases_number.loc[:, 'Total Bases'].apply(integer_to_human_readable,
                                                                  args=(bases_alphabet,)),
                   passed_reads_ratio.apply(ration_to_percentage_string, args=(2,)),
                   passed_bases_number.loc[:, 'Passed Bases'].apply(integer_to_human_readable,
                                                                    args=(bases_alphabet,)),
                   gc_content_ratio.apply(ration_to_percentage_string, args=(2,)),
                   timestamps.loc[:, 'Duration'].apply(lambda x:
                                                        f'{round_to_x_significant(x / 3600, 3)}H')]

    figure = go.Figure(data=[go.Table(header={'values': header_values},
                                      cells={'values': body_values})])

    figure.update_layout(height=100, margin={'r': 5, 'l': 5, 't': 5, 'b': 5})

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


summary_table(csv_files=snakemake.input, output_file=snakemake.output[0])
