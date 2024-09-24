"""Describes function for parsing pore activity file"""

import logging

import pandas as pd

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def parse_pore_activity(input_file, samples, output_files):
    """
    Parses pore activity file and saves to csv file
    :param str input_file: Input file with data to parse
    :param str samples: Samples
    :param list[str] output_files: List of output files names for act, seq and ratio statistics
    :rtype: None
    """

    df = pd.read_csv(input_file)

    df.columns.values[2] = samples

    act_df = (df[df.iloc[:, 0] == "pore"].iloc[:, 2] / 126 / 5000).reset_index(drop=True).T
    seq_df = (df[df.iloc[:, 0] == "strand"].iloc[:, 2] / 126 / 5000).reset_index(drop=True).T

    act_df.replace(0, 1, inplace=True)
    seq_df.replace(0, 1, inplace=True)

    if act_df.shape != seq_df.shape:
        raise ValueError("pore_act and pore_seq have different shapes")

    ratio = (seq_df / act_df) * 100

    act_df.to_csv(output_files[0], header=False)
    seq_df.to_csv(output_files[1], header=False)
    ratio.to_csv(output_files[2], header=False)


try:
    parse_pore_activity(input_file=snakemake.input[0], samples=snakemake.params.samples,
                        output_files=snakemake.output)
except Exception as e:
    logger.exception(e)
    raise e
