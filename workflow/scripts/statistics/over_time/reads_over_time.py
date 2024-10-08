"""Plots number of bases and reads over time"""

import os

import pandas as pd
import plotly.graph_objects as go

from snakemake.script import snakemake

from workflow.scripts.utils import snakemake_file_logger


# pylint: disable=too-many-locals
@snakemake_file_logger
def reads_over_time(path_to_samples, output_file):
    """
    Plots number of reads over time as line plot
    :param list[str] path_to_samples: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    throughput_files = []
    for path_to_sample in path_to_samples:
        for file in os.listdir(path_to_sample):
            if file.startswith('throughput'):
                throughput_files.append(os.path.join(path_to_sample, file))


    figure = go.Figure()

    for path in throughput_files:
        filename = '_'.join(path.split('/')[-4:-2])
        data = pd.read_csv(path).convert_dtypes()
        reads_output_per_minute = data.loc[:, 'Reads'].diff()
        reads_output_per_minute.loc[0] = data.loc[0, 'Reads']
        x_data = data.iloc[:, 0]
        y_data: pd.Series = reads_output_per_minute

        windows_size = y_data.shape[0] // 20
        x_data_averaged = x_data
        y_data_averaged = y_data

        for _ in range(3):
            y_data_averaged = y_data_averaged.rolling(window=windows_size,
                                                      min_periods=windows_size,
                                                      center=True).mean()

        figure.add_trace(go.Scatter(x=x_data,
                                    y=y_data,
                                    name=filename,
                                    legendgroup=filename,
                                    legendgrouptitle={'text': filename}))

        figure.add_trace(go.Scatter(x=x_data_averaged,
                                    y=y_data_averaged,
                                    name=filename + '_averaged',
                                    legendgroup=filename,
                                    legendgrouptitle={'text': filename}))

        figure.update_layout(xaxis_title='time',
                             yaxis_title='number of reads',
                             title='Number of reads over time',
                             legend_groupclick='toggleitem')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False))


reads_over_time(path_to_samples=snakemake.input, output_file=snakemake.output[0])
