"""Creates mapping table summary"""


import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_records,
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
    header_values = ['Experiment', 'Sequencing Efficiency']
    body_values = [[], []]
    for file in sam_files:
        counter = 0
        with open(file, encoding='utf8') as f:
            for _ in range(4):
                f.readline()

            group_of_reads = []
            end_of_file = False
            while True:
                if end_of_file:
                    break
                while True:
                    read = f.readline().strip()
                    if not read:
                        end_of_file = True
                        break
                    if not group_of_reads:
                        group_of_reads.append(read.split())
                    if read.split()[0] == group_of_reads[0][0]:
                        group_of_reads.append(read.split())
                    else:
                        break

                if len(group_of_reads) == 2:
                    read1, read2 = parse_sam_records(group_of_reads)
                    if read1[11]['NM'] <= 3 and read2[11]['NM'] <= 3 and read1[2] != read2[2]:
                        is_reverse1 = bin(int(read1[1]))[2:].rjust(5, '0')[-5]
                        is_reverse2 = bin(int(read2[1]))[2:].rjust(5, '0')[-5]
                        if is_reverse1 == is_reverse2:
                            counter += 1
                group_of_reads = [read.split()]

        body_values[0].append(file.split('/')[-1].split('.')[0])
        body_values[1].append(integer_to_human_readable(counter))

    figure = go.Figure(data=[go.Table(header={'values': header_values},
                                      cells={'values': body_values})])

    figure.update_layout(height=100, margin={'l': 5, 'r': 5, 't': 5, 'b': 5})


    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


mapping_table(sam_files=snakemake.input, output_file=snakemake.output[0])
