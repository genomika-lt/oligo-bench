"""Describes function for parsing samtools statistics"""

import logging

import pandas as pd

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def parse_samtools_stats(input_file, samples, output_file):
    """
    Parses samtools statistics input file and save results to output_file.
    :param str input_file: Input file with data to parse
    :param str samples: Samples
    :param list[str] output_file: List of output file for mq, in, dl and co statistics
    :rtype: None
    """

    mq_lines = []
    id_lines = []
    co_lines = []

    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            values = line.strip().split('\t')
            if line.startswith("MAPQ"):
                mq_lines.append(values)
            elif line.startswith("ID"):
                id_lines.append(values)
            elif line.startswith("COV"):
                co_lines.append(values)
            else:
                raise ValueError("Unexpected line in samtools statistics file")

    mq_df = pd.DataFrame(mq_lines, columns=["Type", "MAPQ", samples]).iloc[:, [2]].T
    in_df = pd.DataFrame(id_lines, columns=["Type", "Length", samples, "Num_Del"]).iloc[:, [2]].T
    dl_df = pd.DataFrame(id_lines, columns=["Type", "Length", "Num_Ins", samples]).iloc[:, [3]].T
    co_df = pd.DataFrame(co_lines, columns=["Type", "Range", "Depth", samples]).iloc[:, [3]].T

    mq_df.to_csv(output_file[0], header=False)
    in_df.to_csv(output_file[1], header=False)
    dl_df.to_csv(output_file[2], header=False)
    co_df.to_csv(output_file[3], header=False)


try:
    parse_samtools_stats(input_file=snakemake.input[0],
                         samples=snakemake.params.samples,
                         output_file=snakemake.output)
except Exception as e:
    logger.exception(e)
    raise e
