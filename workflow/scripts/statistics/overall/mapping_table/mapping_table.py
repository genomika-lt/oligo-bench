"""Plots summary table for samples from precalculated csv files and saves to html. """


import pandas as pd
import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (snakemake_file_logger,
                           round_to_x_significant,
                           integer_to_human_readable)

def ratio_to_percentage_string(ratio: float, digits_after_point=2) -> str:
    """
    Converts ratio to percantage string.
    :param ratio: ratio to convert.
    :param digits_after_point: number of digits after the point.
    :return: String with percentage and % at the end.
    """
    return f'{100 * round_to_x_significant(ratio, 2 + digits_after_point):.{digits_after_point}f}%'

# pylint: disable=too-many-locals
@snakemake_file_logger
def mapping_table(csv_files, output_file):
    """
    Plots mapping table for samples from precalculated csv files and saves to html.
    :param dict[str, str] csv_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    bases_alphabet = ('b', 'Kb', 'Mb', 'Gb', 'Tb', 'Pb')

    forward_reads_ratio = pd.read_csv(csv_files['forward_reads'])
    total_reads = pd.read_csv(csv_files['total_reads'])
    mapped_reads = pd.read_csv(csv_files['mapped_reads'])
    total_bases = pd.read_csv(csv_files['total_bases'])
    mapped_bases = pd.read_csv(csv_files['mapped_bases'])

    mapped_bases_percent = mapped_bases.loc[:, 'Mapped Bases'] / total_bases.loc[:, 'Total Bases']
    mapped_reads_percent = mapped_reads.loc[:, 'Mapped Reads'] /  total_reads.loc[:, 'Total Reads']

    header_values = ['Sample ID',
                     'Mapped Reads',
                     'Mapped Bases',
                     'Forward Reads']

    body_values_number = [forward_reads_ratio.loc[:, 'Sample'],
                   mapped_reads.loc[:, 'Mapped Reads'].apply(integer_to_human_readable),
                   mapped_bases.loc[:, 'Mapped Bases'].apply(integer_to_human_readable,
                                                             args=(bases_alphabet,)),
                   forward_reads_ratio.loc[:, 'Forward Reads'].apply(ratio_to_percentage_string, args=(2,))]

    body_values_percentage = [forward_reads_ratio.loc[:, 'Sample'],
                    mapped_reads_percent.apply(ratio_to_percentage_string, args=(3,)),
                    mapped_bases_percent.apply(ratio_to_percentage_string, args=(3,)),
                    forward_reads_ratio.loc[:, 'Forward Reads'].apply(ratio_to_percentage_string, args=(2,))]

    figure = go.Figure(data=[go.Table(header={'values': header_values},
                                      cells={'values': body_values_number,
                                             'height': 25})])
    figure.update_layout(height=25 * forward_reads_ratio.shape[0] + 50,
                         margin={'r': 5, 'l': 5, 't': 5, 'b': 5})

    figure.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                active=1,
                buttons=list([
                    dict(label="Percentage",
                         method="update",
                         args=[{'cells': {'values': body_values_percentage}}]),
                    dict(label="Quantity",
                         method="update",
                         args=[{'cells': {'values': body_values_number}}]),
                ]),
            )
        ])

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


mapping_table(csv_files=snakemake.input, output_file=snakemake.output[0])
