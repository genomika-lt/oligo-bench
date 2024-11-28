"""Creates mapping table summary"""


import plotly.graph_objects as go
import pysam

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (group_sam_bam_file,
                           snakemake_file_logger,
                           integer_to_human_readable)


@snakemake_file_logger
def filter_2_primers(sam_file, output_file):
    """
    Plots number of bases and reads over time as line plot and cumulative output
    :param list[str] sam_files: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """


    allowed_primer_errors = snakemake.config['allowed_primer_errors']

    # pylint: disable=no-member
    read_file = pysam.AlignmentFile(sam_file, 'r', check_sq=False)
    # pylint: disable=no-member
    write_file = pysam.AlignmentFile(output_file, "w", template=read_file, check_sq=False)
    for group in group_sam_bam_file(sam_file):
        if len(group) == 2:
            read1, read2 = group

            are_reverse_and_forward = read1.reference_name != read2.reference_name
            are_mapped_correctly = read1.is_forward == read2.is_forward

            if read1.reference_name != 'Forward5':
                read1, read2 = read2, read1
            is_number_of_errors_allowed = (read1.get_tag('NM') <= allowed_primer_errors
                          and read2.get_tag('NM') <= allowed_primer_errors)

            if are_reverse_and_forward and are_mapped_correctly and is_number_of_errors_allowed:
                write_file.write(read1)
                write_file.write(read2)


filter_2_primers(sam_file=snakemake.input[0], output_file=snakemake.output[0])
