"""Describes function for counting mismatches in bam file"""


import logging

import pysam
import pandas as pd

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def count_mism_bam(bam_file, output_file, threads):
    """
    Counts mismatches in bam and writes output to out
    :param str bam_file: Input bam file
    :param str output_file: Output file to save data
    :param int threads: Number of threads to use
    :rtype: None
    """

    pysam.index(bam_file)

    inp = pysam.AlignmentFile(bam_file, 'rb', threads=threads)

    # init a dict to count mismatches
    mism_counts = {
        'NM0': 0, 'NM1': 0, 'NM2': 0,
        'NM3': 0, 'NM4': 0, '>NM5': 0
    }

    # parse the BAM file
    for read in inp.fetch():
        # use primary reads
        if not read.is_supplementary and not read.is_secondary:
            number_of_mismatches = read.get_tag("NM")
            if number_of_mismatches < 5:
                mism_counts['NM' + str(number_of_mismatches)] += 1
            else:
                mism_counts['>NM5'] += 1

    inp.close()

    # flatten
    mism_df = pd.DataFrame(list(mism_counts.items()),
                           columns=['mism', 'count'])

    # write output
    mism_df.to_csv(output_file, header=False, index=False)


try:
    count_mism_bam(bam_file=snakemake.input[0], output_file=snakemake.output[0],
                   threads=snakemake.threads)
except Exception as e:
    logger.exception(e)
    raise e
