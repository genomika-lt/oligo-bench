"""Counts total number of reads"""


from snakemake.script import snakemake

from workflow.scripts.utils import file_logger


@file_logger
def finalise_report(reports, output_file):
    """
    Counts total number of reads across samples
    :param str reports: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    with open(output_file, 'w', encoding='utf-8') as g:
        for report in reports:
            with open(report, 'r', encoding='utf-8') as f:
                g.write(f.read())


finalise_report(reports=snakemake.input, output_file=snakemake.output[0])
