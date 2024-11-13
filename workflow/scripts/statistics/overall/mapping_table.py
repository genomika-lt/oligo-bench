"""Creates mapping table summary"""


import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (group_sam_bam_file,
                           snakemake_file_logger,
                           integer_to_human_readable)


@snakemake_file_logger
def mapping_table(sam_files, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] sam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    # pylint: disable=too-many-locals
    header_values = ['Experiment', 'Mapped 2 primers']
    body_values = [[], []]
    allowed_primer_errors = snakemake.config['allowed_primer_errors']
    for file in sam_files:
        counter = 0
        for group in group_sam_bam_file(file):
            if len(group) == 2:
                read1, read2 = group

                condition1 = read1.reference_name != read2.reference_name
                condition2 = read1.is_forward == read2.is_forward
                condition3 = (read1.get_tag('NM') <= allowed_primer_errors
                              and read2.get_tag('NM') <= allowed_primer_errors)

                if condition1 and condition2 and condition3:
                    counter += 1

        body_values[0].append(file.split('/')[-1].split('.')[0])
        body_values[1].append(integer_to_human_readable(counter))

    figure = go.Figure(data=[go.Table(header={'values': header_values},
                                      cells={'values': body_values})])

    figure.update_layout(height=100, margin={'l': 5, 'r': 5, 't': 5, 'b': 5})


    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


mapping_table(sam_files=snakemake.input, output_file=snakemake.output[0])
