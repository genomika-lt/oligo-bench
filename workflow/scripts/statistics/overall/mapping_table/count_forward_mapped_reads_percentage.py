"""Creates mapping table summary"""


import plotly.graph_objects as go
import pysam
from pandas import DataFrame
from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_bam_file,
                                snakemake_file_logger,
                                integer_to_human_readable)


@snakemake_file_logger
def count_forward_mapped_reads_percentage(sam_bam_files, output_file):
    """
    Counts forward and reversed mapped reads from sam files and gives ratio.
    :param list[str] sam_bam_files: List with paths to SAM/BAM files.
    :param str output_file: Output file to save data.
    :rtype: None
    """


    data = {'Sample': [],
            'Forward Reads': []}

    for file in sam_bam_files:
        forward_counter = 0
        reverse_counter = 0

        for read in parse_sam_bam_file(file):
            if read.is_forward:
                forward_counter += 1
            else:
                reverse_counter += 1

        sample = file.split('/')[-1].split('.')[0][16:]
        data['Sample'].append(sample)
        data['Forward Reads'].append(forward_counter / (forward_counter + reverse_counter))

    data_frame = DataFrame(data)
    data_frame.to_csv(output_file, index=False)


count_forward_mapped_reads_percentage(sam_bam_files=snakemake.input, output_file=snakemake.output[0])
