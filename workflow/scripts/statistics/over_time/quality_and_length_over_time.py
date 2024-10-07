"""Plots dependency between quality and length of reads"""


from datetime import datetime

import pysam
import pandas as pd
import plotly.express as px

from snakemake.script import snakemake

from workflow.scripts.utils import snakemake_file_logger, parse_sam_records


@snakemake_file_logger
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
        parsed_records = parse_sam_records(records)

        filename = path_to_sample.split('/')[-1][:-4]

        parsed_records.sort(key=lambda x: x[11]['st'])

        first_read_time = datetime.strptime(parsed_records[0][-1]['st'], '%Y-%m-%dT%H:%M:%S.%f+00:00')

        length_quality_time_data = {'Quality': [],
                                    'Length': [],
                                    'Time': [],
                                    'Sample': []}

        for record in parsed_records:
            length_quality_time_data['Quality'].append(record[-1]['qs'])
            length_quality_time_data['Length'].append(len(record[9]))
            current_read_time = datetime.strptime(record[-1]['st'],'%Y-%m-%dT%H:%M:%S.%f+00:00')
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


quality_and_length_over_time(bam_files=snakemake.input, output_file=snakemake.output[0])
