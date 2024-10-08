"""Plots mean reads length over time as line"""


from datetime import datetime, timedelta

import pysam
import pandas as pd
import plotly.express as px

from snakemake.script import snakemake

from workflow.scripts.utils import parse_sam_records, snakemake_file_logger


# pylint: disable=too-many-locals
@snakemake_file_logger
def read_length_over_time_as_line(bam_files, output_file):
    """
    Plots gc distribution over time in samples and saves to html
    :param list[str] bam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    plotting_data = pd.DataFrame()

    for file in bam_files:
        out = pysam.view('-o', 'out.sam', file)
        records = [record.split() for record in  out.split('\n')[:-1]]
        parsed_records = parse_sam_records(records)
        parsed_records.sort(key=lambda x: x[11]['st'])

        read_length_over_minute = []
        minutes_counter = 1
        bases_counter = 0
        number_of_reads_counter = 0
        first_read_time = datetime.strptime(parsed_records[0][11]['st'],
                                            '%Y-%m-%dT%H:%M:%S.%f+00:00')

        for record in parsed_records:
            current_read_time = datetime.strptime(record[11]['st'], '%Y-%m-%dT%H:%M:%S.%f+00:00')
            read = record[9]
            if current_read_time - first_read_time < timedelta(minutes=minutes_counter):
                bases_counter += len(read)
                number_of_reads_counter += 1
            else:
                read_length_over_minute.append(bases_counter / number_of_reads_counter)
                minutes_counter += 1
                number_of_reads_counter = 1
                bases_counter = len(read)

        read_length_over_minute.append(bases_counter / number_of_reads_counter)

        plotting_data = plotting_data.join(pd.DataFrame(read_length_over_minute,
                                                        columns=[file.split('/')[-1][:-4]]),
                                           how='outer')

    figure = px.line(plotting_data,
                     x=plotting_data.index,
                     y=plotting_data.columns,
                     title='Read length during time'
                     )

    figure.update_layout(xaxis_title='Time (Minutes)',
                         yaxis_title='Mean Length Of Read',
                         legend_title='Samples')

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


read_length_over_time_as_line(bam_files=snakemake.input, output_file=snakemake.output[0])
