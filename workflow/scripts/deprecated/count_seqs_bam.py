"""Describes function for counting total number of sequences in bam file"""


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


def count_seqs_bam(bam, ref, samples, output_file, threads):
    """
    Counts total number of sequences in bam file
    :param str bam: Input bam file
    :param str ref: Reference file
    :param str samples: Samples
    :param str output_file:  Output file to save data
    :param int threads: Number of threads to use
    :rtype: None
    """

    pysam.index(bam)

    inp = pysam.AlignmentFile(bam, 'rb', threads=threads)

    # get reference names
    fa = pysam.FastaFile(ref)
    fa_names = fa.references

    # init a dict to count sequences
    ref_counts = {name: 0 for name in fa_names}

    # parse bam
    for read in inp.fetch():
        # use primary reads
        if not read.is_supplementary and not read.is_secondary:
            ref_name = read.reference_name
            if ref_name in ref_counts:
                ref_counts[ref_name] += 1

    inp.close()

    # flatten
    ds = pd.Series(ref_counts)

    # calculate freq
    freq = ds.value_counts().sort_index()

    # fill in zeros
    rang = range(ds.min(), ds.max() + 1)
    freq = freq.reindex(rang, fill_value=0)

    # convert
    freq_df = pd.DataFrame([freq], index=[samples])

    # write output
    freq_df.to_csv(output_file, header=True)

try:
    count_seqs_bam(bam=snakemake.input[0], ref=snakemake.input[1], samples=snakemake.params.samples,
                   output_file=snakemake.output[0], threads=snakemake.threads)
except Exception as e:
    logger.exception(e)
    raise e
