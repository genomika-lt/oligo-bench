"""Plots pore activity over pores time for each sample"""

import logging
import os

import pandas as pd
import plotly.graph_objects as go

from workflow.scripts.utils import file_logger

from snakemake.script import snakemake


@file_logger
def pore_activity(path_to_samples, output_file):
    """
    Plots pore activity based on minknow output csv file
    :param list[str] path_to_samples: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    pore_activity_files = []
    for path in path_to_samples:
        for file in os.listdir(path):
            if file.startswith("pore_activity"):
                pore_activity_files.append(os.path.join(path, file))

    figure = go.Figure()
    counter = 0
    for file in pore_activity_files:
        counter += 1
        filename = '_'.join(file.split('/')[-4:-2])

        data = pd.read_csv(file, dtype={'Channel State': 'str',
                                          'Experiment Time (minutes)': 'int',
                                          'State Time (samples)': 'float'})
        data.loc[:, 'activity'] = data.iloc[:, 2] / 126 / 5000

        meta_data = {'rank': {'strand': 1,
                              'adapter': 2,
                              'pore': 3,
                              'no_pore': 4},
                     'color': {'strand': '#00ff00',
                               'adapter': '#ede797',
                               'pore': '#00cc00',
                               'no_pore': '#4da9c3'},
                     'visible': {'strand': True,
                                 'adapter': True,
                                 'pore': True,
                                 'no_pore': 'legendonly'}
                     }

        groups = data.groupby('Channel State')
        for group in groups.groups:
            figure.add_trace(go.Scatter(x=groups.get_group(group)['Experiment Time (minutes)'],
                                        y=groups.get_group(group)['activity'],
                                        name=filename,
                                        legendgroup=group,
                                        legendgrouptitle={'text': group},
                                        legendrank=meta_data['rank'].get(str(group), 1000),
                                        visible=meta_data['visible'].get(str(group), 'legendonly'),
                                        marker={'color': meta_data['color'].get(str(group), '#000000')}))


    figure.update_layout(legend_groupclick='toggleitem')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


pore_activity(path_to_samples=snakemake.input, output_file=snakemake.output[0])
