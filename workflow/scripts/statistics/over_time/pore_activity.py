"""Plots pore activity over pores time for each sample"""

import os
import json

import pandas as pd
import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger



@snakemake_file_logger
def pore_activity(path_to_samples, output_file):
    """
    Plots pore activity based on minknow output csv file
    :param list[str] path_to_samples: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    files = []
    for path in path_to_samples:
        csv_file = ''
        json_file = ''
        for file in os.listdir(path):
            if file.startswith("pore_activity"):
                csv_file = os.path.join(path, file)
            if file.endswith(".json"):
                json_file = os.path.join(path, file)
        files.append((csv_file, json_file))

    figure = go.Figure()

    for csv_file, json_file in files:
        sample_name = csv_file.split('/')[-3]

        with (open(json_file, 'r', encoding='utf-8') as f):
            json_data = json.load(f)
            sample_rate = (json_data['acquisitions'][0]['acquisition_run_info']
                                       ['config_summary']['sample_rate'])
            number_of_wells = (json_data['protocol_run_info']['flow_cell']['channel_count'] *
                               json_data['protocol_run_info']['flow_cell']['wells_per_channel'])

        data = pd.read_csv(csv_file, dtype={'Channel State': 'str',
                                          'Experiment Time (minutes)': 'int',
                                          'State Time (samples)': 'float'})

        data.loc[:, 'activity'] = (data.loc[:, 'State Time (samples)'] / sample_rate / 60 /
                                   number_of_wells)

        groups = data.groupby('Channel State')
        figure.add_trace(go.Scatter(x=groups.get_group('strand')['Experiment Time (minutes)'],
                                    y=groups.get_group('strand')['activity'],
                                    name=sample_name))

    figure.update_layout(xaxis_title='Time in minutes',
                         yaxis_title='Number of pores',
                         title='Sequencing Pores Activity',
                         legend_title='Samples')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


pore_activity(path_to_samples=snakemake.input, output_file=snakemake.output[0])
