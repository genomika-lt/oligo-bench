"""Plots summary table and saves to html"""


import pysam

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger


@snakemake_file_logger
def filter_basecalled(bam_file, output_file):
    """
    Plots gc distribution over time in samples and saves to html
    :param str bam_file: Folder with fastq reads
    :param str output_file: Output file to save data to
    :rtype: None
    """

    # pylint: disable=no-member
    data = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
    # pylint: disable=no-member
    filtered = pysam.AlignmentFile(output_file, "wb", template=data, check_sq=False)
    for read in data:
        for tag in read.tags:
            if tag[0] == 'qs' and tag[1] > 9:
                filtered.write(read)


filter_basecalled(bam_file=snakemake.input[0], output_file=snakemake.output[0])
