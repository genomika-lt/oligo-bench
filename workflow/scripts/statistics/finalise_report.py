"""Counts total number of reads"""

import logging

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


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


try:
    finalise_report(reports=snakemake.input, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
