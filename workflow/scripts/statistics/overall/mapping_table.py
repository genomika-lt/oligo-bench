"""Creates mapping table summary"""


import pysam
import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from ....scripts.utils import parse_sam_records, snakemake_file_logger



@snakemake_file_logger
def mapping_table(sam_files, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] sam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    #
    # for file in sam_files:
    #     for read in parse_sam_records(file):
    #         pass
    # with open(output_file, 'w', encoding='utf-8') as g:
    #     g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


mapping_table(sam_files=snakemake.input, output_file=snakemake.output[0])
