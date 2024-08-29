"""Plots pore activity based on minknow output csv file"""

import logging

import pandas as pd
import plotly.express as px

from snakemake.script import snakemake

# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def pore_activity(pore_activity_files, output_file):
    """
    Plots pore activity based on minknow output csv file
    :param list[str] pore_activity_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    plotting_data = pd.DataFrame()

    for file in pore_activity_files:
        groups = pd.read_csv(file).groupby('Channel State')

        sequencing = groups.get_group('strand')
        sequencing.reset_index(drop=True, inplace=True)

        sequencing['activity'] = pd.to_numeric(sequencing.iloc[:, 2], downcast='float') / 500 / 60
        sequencing['sample'] = '_'.join(file.split('/')[-3:-1]) + '_sequencing'

        plotting_data = pd.concat([plotting_data, pd.DataFrame(sequencing)], axis=0)

    figure = px.line(plotting_data,
                     x='Experiment Time (minutes)',
                     y='activity',
                     title='pore activity during time',
                     color='sample')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False))

try:
    pore_activity(pore_activity_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
