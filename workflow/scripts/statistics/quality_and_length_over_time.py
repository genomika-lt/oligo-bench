"""Plots dependency between quality and length of reads"""

import logging
from datetime import datetime

import pysam
import pandas as pd
import plotly.express as px

from snakemake.script import snakemake

# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def quality_and_length_over_time(bam_files, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    plotting_data = pd.DataFrame()

    for path_to_sample in bam_files:
        out = pysam.view('-o', 'out.sam', path_to_sample)
        records = [line.split() for line in out.split('\n')[:-1]]

        filename = path_to_sample.split('/')[-1]
        date_index = quality_index = -1
        read_index = 9
        for i in range(len(records[0])):
            if records[0][i].startswith('st:Z:'):
                date_index = i
            if records[0][i].startswith('qs:i:'):
                quality_index = i

        records.sort(key=lambda x: x[date_index])

        first_read_time = datetime.strptime(records[0][date_index], 'st:Z:%Y-%m-%dT%H:%M:%S.%f+00:00')

        length_quality_time_data = {'Quality': [],
                                    'Length': [],
                                    'Time': [],
                                    'Sample': []}

        for record in records:
            length_quality_time_data['Quality'].append(int(record[quality_index][5:]))
            length_quality_time_data['Length'].append(len(record[read_index]))
            current_read_time = datetime.strptime(record[date_index],'st:Z:%Y-%m-%dT%H:%M:%S.%f+00:00')
            length_quality_time_data['Time'].append((current_read_time - first_read_time).seconds // 60)
            length_quality_time_data['Sample'].append(filename)

        plotting_data = pd.concat([plotting_data, pd.DataFrame(length_quality_time_data).convert_dtypes()])

    figure = px.density_contour(plotting_data,
                                x='Length',
                                y='Quality',
                                animation_frame='Time',
                                marginal_x='histogram',
                                marginal_y='histogram',
                                color='Sample')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))

try:
    quality_and_length_over_time(bam_files=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
