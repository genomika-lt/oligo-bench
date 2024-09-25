"""Plots gc distribution over time in samples and saves to html"""


from datetime import datetime, timedelta

import pysam
import pandas as pd
import plotly.express as px

from snakemake.script import snakemake

from workflow.scripts.utils import file_logger


@file_logger
def gc_over_time(bam_files, output_file):
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
        records.sort(key=lambda x: x[17])

        # Counting GC% over time
        gc_distributions_over_minute = []
        gc_counter = 0
        minutes_counter = 1
        total_bases_counter = 0
        first_read_time = datetime.strptime(records[0][17], 'st:Z:%Y-%m-%dT%H:%M:%S.%f+00:00')

        for record in records:
            current_read_time = datetime.strptime(record[17], 'st:Z:%Y-%m-%dT%H:%M:%S.%f+00:00')
            read = record[9]
            if current_read_time - first_read_time < timedelta(minutes=minutes_counter):
                gc_counter += read.count('G') + read.count('C')
                total_bases_counter += len(read)
            else:
                gc_distributions_over_minute.append(gc_counter / total_bases_counter)
                minutes_counter += 1
                gc_counter = read.count('G') + read.count('C')
                total_bases_counter = len(read)
        gc_distributions_over_minute.append(gc_counter / total_bases_counter)

        plotting_data = plotting_data.join(pd.DataFrame(gc_distributions_over_minute, columns=[file.split('/')[-1]]), how='outer')

    plotting_data['Time (Minutes)'] = plotting_data.index

    figure = px.line(plotting_data, x='Time (Minutes)', y=plotting_data.columns, title='GC/bases during time')
    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


gc_over_time(bam_files=snakemake.input, output_file=snakemake.output[0])
